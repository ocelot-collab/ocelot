//#include "gauss_cell_integr.h"
#include <omp.h>
#include <math.h>
//template <typename T> inline constexpr // for gcc
template <typename T> inline
int sign(T val) {
    return (T(0) < val) - (val < T(0));
}
inline double erf_my(double x)
{  
	/* erf(z) = 2/sqrt(pi) * Integral(0..x) exp( -t^2) dt
	erf(0.01) = 0.0112834772 erf(3.7) = 0.9999998325
	Abramowitz/Stegun: p299, |erf(z)-erf| <= 1.5*10^(-7) 
	*/
	//int signx = sign(x);
	double y = 1.0 / ( 1.0 + 0.3275911 * fabs(x));   
	return sign(x)*(1. - (((((
        + 1.061405429  * y
        - 1.453152027) * y
        + 1.421413741) * y
        - 0.284496736) * y 
        + 0.254829592) * y)*exp (-x * x));      
}

int data_extention_fast(const double *data, const int nx, const int ny, double *data_ext)
{
    //# vertical -> y -> i
    //# horizontal -> x -> j    

    int nx_e = 2 + nx;
    int ny_e = 2 + ny;

    for (int i = 0; i< ny; i++)
        for (int j = 0 ; j< nx; j++)
            data_ext[(i+1)*nx_e + j+1] = data[i*nx +j];
    // top and bottom
    for (int j = 0 ; j< nx; j++)
    {
        data_ext[j+1] = data[j];
        data_ext[(ny_e-1)*nx_e + j+1] = data[(ny-1)*nx + j];
    }
    // left and right
    for (int i = 0; i < ny_e ; i++)
    {
        data_ext[i*nx_e ] = data_ext[i*nx_e +1];
        data_ext[i*nx_e + nx_e-1] = data_ext[i*nx_e +nx_e-2];
    }
    return 1;
}

int data_extention_1D(const double *data, const int nx, double *data_ext)
{
    //# vertical -> y -> i
    //# horizontal -> x -> j

    int nx_e = 2 + nx;
    data_ext[0] = data[0];
    data_ext[nx_e - 1] = data[nx-1];
    for (int j = 0 ; j< nx; j++)
        data_ext[j+1] = data[j];
    return 1;
}

inline void precalculate_XY2(const int nx, const double *X, const double sx, const double Tlim,
	double *Xnew, double *expX, double *erfX )
{
	int nx_e = nx + 2;
    double sq2 = sqrt(2.);
    double sq2sx = sq2*sx;
	
    double *Xcalc = new double [nx+nx - 1];
    double *expXcalc = new double [nx+nx - 1];
    double *erfXcalc = new double [nx+nx - 1];
    Xnew[nx_e - 1] = Tlim;
    double argx0 = Xnew[nx_e - 1]/sq2sx;
    expX[nx_e-1] = exp(-argx0*argx0);
    erfX[nx_e-1] = erf_my(argx0);
    //#expX[0,-1] = 0.
    //#erfX[0, -1] = 1
    int p = 0;
    for (int n = 0; n<2; n++)// in range(2):
	{
        for (int j = n; j<nx; j++)// in range(n,nx):
		{
            //Xcalc[p] = X[j] - X[(nx-1)*n];
            //double argx0 = Xcalc[p]/sq2sx;
            //expXcalc[p] = exp(-argx0*argx0);
            //erfXcalc[p] = erf(argx0);
            //# extra function
            Xnew[j*nx_e+(nx_e-1)*n] = (-1 + 2*n)*Tlim;
            expX[j*nx_e+(nx_e-1)*n] = expX[nx_e-1];
            erfX[j*nx_e+(nx_e-1)*n] = (-1. + 2.*n)*erfX[nx_e-1];
            p = p +1;
		}
	}

    for(int i = 0; i<nx;i++)
	{
        for(int j = 0; j<nx; j++)
		{
            Xnew[i*nx_e+j+1] = X[j] - X[i];
            double argx0 = Xnew[i*nx_e+j+1] /sq2sx;
            expX[i*nx_e+j+1] = exp(-argx0*argx0);
            erfX[i*nx_e+j+1] = erf_my(argx0);
		}
	}
	delete [] Xcalc;
	delete []expXcalc;
	delete [] erfXcalc;
    //return Xnew,  expX,  erfX
}

inline double EMITxDATA_fast2(const double *data_ext, const int nx, const int ny,
	const double *Xnew, const double *Ynew, const double sx, const double sy,
	const double *expX, const double *erfX, const double *expY, const double *erfY)
{
    //# vertical -> y -> i
    //# horizontal -> x -> j
	double pi = 3.14159265358979324;
	double sq_2_pi = sqrt(2.*pi);
    double integral = 0.;
	double pi4 = 4.*pi;
	double ky = sq_2_pi*sy;
	double kx = sq_2_pi*sx;
	double kk = 2.*sx*sy; 
	// initialization 
	//double f00, f10, f01,f11, x0, x1, a,b,c,d,A,B,d_expX,d_erfX, den, dF1, dF2;
    for( int i = 0; i < (ny-1); i++) // for exp(Y)
    {    
        double y0 = Ynew[i]; 
        double y1 = Ynew[i + 1];

        double d_erfY = erfY[i] - erfY[i+1];
        double d_expY = expY[i] - expY[i+1];

		double kx_d_erfY = kx*d_erfY ;
		double d_expY_kk = d_expY*kk;
		double d_expY_ky = d_expY*ky;
		double pi_d_erfY = pi*d_erfY;
		//int i_nx = i*nx;
		for(int j = 0; j<nx-1; j++)
		{
       
			// put integration over cell
			double f00 = data_ext[i*nx+j];
			double f10 = data_ext[i*nx+j+1]; 
			double f01 = data_ext[i*nx+j+nx];
			double f11 = data_ext[i*nx+j + nx+1];

			double x0 = Xnew[j];
			double x1 = Xnew[j + 1];

			double d_expX = expX[j+1] - expX[j];
			double d_erfX = erfX[j+1] - erfX[j];
			double den = (x0 - x1)*(y0 - y1);
			double a = (f11*x0 - f01*x1)*y0 - (f10*x0 - f00*x1)*y1;
			double dF1 = f00 - f10;
			double dF2 = f01 - f11;
			double b = dF2*y0 - dF1*y1;
			double c = (f10 - f11)*x0 - (f00 - f01)*x1;
			double d = dF1 - dF2;

			double A = b*kx_d_erfY - d*d_expY_kk;
			double B = c*d_expY_ky - a*pi_d_erfY;

			integral += (A*d_expX + B*d_erfX)/(pi4*den);
			//put integration over cell 
		}
    }
    return integral;
}

inline double Gauss_x_Func_1D(const double *data_ext,const int nx_e, const double sx, const double *Xnew, const double *expX, const double *erfX)
{
    double pi = 3.14159265358979324;
    double sq_2_pi = sqrt(2.*pi);
    double integral = 0.;
    for (int n = 0; n<nx_e-1; n++)
    {

        double X0 = Xnew[n];
        double X1 = Xnew[n+1];

        double f0 = data_ext[n];
        double f1 = data_ext[n+1];

        double den = X0 - X1;
        double a = f1*X0 - f0*X1;
        double b = f0 - f1;
        double erfX0 = erfX[n];
        double erfX1 = erfX[n+1];

        double expX0 = expX[n];
        double expX1 = expX[n+1];

        double A = b*sx/sq_2_pi*(expX0 - expX1);
        double B = a/2.*(-erfX0+erfX1);
        integral += (A+B)/den;
    }
    return integral;
}

extern "C"
{
int conv_1D(double *screen, const double *X,  const int nx, const int nx_add, const double sx)
{
    double *data_ext = new double [nx + 2];
    data_extention_1D(screen, nx, data_ext);
    int nx_e = nx+2;
    int ret = 0;
    double *Xnew = new double [nx_e*nx];
    double *expX = new double [nx_e*nx];
    double *erfX = new double [nx_e*nx];
    double Xlim = fabs(X[nx-1] - X[0]);
    double Tlim = Xlim*1000.;

    precalculate_XY2(nx, X, sx, Tlim, Xnew, expX, erfX );

    for (int j = nx_add; j < nx - nx_add; j++)
    {
        screen[j] = Gauss_x_Func_1D(data_ext, nx_e, sx, &Xnew[nx_e*j], &expX[nx_e*j], &erfX[j*nx_e]);
    }
    delete [] Xnew;
    delete [] expX;
    delete [] erfX;
    delete [] data_ext;
    return ret;
}
int conv_2D(double *screen, const double *X, const double *Y, const int nx, const int ny, const int nx_add, const int ny_add, const double sx, const double sy)
{
    //# vertical -> y -> i
    //# horizontal -> x -> j
    int ret = 0;
    int nx_e = nx + 2;
    int ny_e = ny + 2;
    double *data_ext = new double [ny_e*nx_e];
    data_extention_fast(screen, nx, ny,data_ext);
    double Xlim = fabs(X[nx-1] - X[0]);
    double Ylim = fabs(Y[ny-1] - Y[0]);
    double Tlim = 1000.;
    if (Xlim>Ylim)
        Tlim = Xlim*1000.;
    else
        Tlim = Ylim*1000.;
    double *Ynew = new double [ny_e*ny];
    double *expY = new double [ny_e*ny];
    double *erfY = new double [ny_e*ny];
    double *Xnew = new double [nx_e*nx];
    double *expX = new double [nx_e*nx];
    double *erfX = new double [nx_e*nx];

    precalculate_XY2( nx,   X,  sx, Tlim,  Xnew, expX,erfX  );
    precalculate_XY2( ny,   Y,  sy, Tlim,  Ynew, expY,erfY  );
    #pragma omp parallel for
    for (int i = ny_add;  i < ny - ny_add; i++)
    {
        for (int j = nx_add; j < nx - nx_add; j++)
        {
            screen[i*nx+j]  = EMITxDATA_fast2(data_ext, nx_e, ny_e, &Xnew[nx_e*j], &Ynew[ny_e*i],  sx,  sy,
                &expX[nx_e*j], &erfX[j*nx_e], &expY[ny_e*i], &erfY[ny_e*i]);
        }
    }
    delete [] Xnew;
    delete [] expX;
    delete [] erfX;
    delete [] Ynew;
    delete [] expY;
    delete [] erfY;
    delete [] data_ext;
    return ret; //erf(-23)*1000;
}


}



