#include "spline.h"
//#include <vector>

//using namespace std;
void derivat(const double *arZ,
             const double *arFz,
             const int size,
             double *deriv1,
             double *deriv2)
{
    int N=size-1;
    double h0;
    double h1;
    double h2;
    double df;
    double F0;
    double F1;
    double a11;
    double a12;
    double a21;
    double a22;
    double det;
    double a0;
    double b0;
    double c0;
    double derivative = 0.;
    int n = 0;
    int i = 0;
    for(i=0; i<2;i++)
    {
        n = i*(N-3);
        h0 = arZ[n+1] - arZ[n];
        h1 = arZ[n+2] - arZ[n];
        h2 = arZ[n+3] - arZ[n];
        df = (arFz[n + 1] - arFz[n])/h0;
        F0 = (arFz[n + 2] - arFz[n])/h1 - df;
        F1 = (arFz[n + 3] - arFz[n])/h2 - df;
        a11 = h1*h1 - h0*h0;
        a12 = h1 - h0;
        a21 = h2*h2 - h0*h0;
        a22 = h2 - h0;
        det = a11*a22 - a12*a21;
        a0 = (a22*F0 - a12*F1)/det;
        b0 = (-a21*F0 + a11*F1)/det;
        c0 = df - a0*h0*h0 - b0*h0;

        derivative = (3.*a0*(arZ[N*i]-arZ[n*i])*(arZ[N*i]-arZ[n*i]) + 2.*b0*(arZ[N*i]-arZ[n*i]))*i + c0;
        if(i==0)
        {
            *deriv1 = derivative;
        }
        else
        {
            *deriv2 = derivative;
        }

    }
}



void moment( const double *arZ,
             const double *arFz,
             const int size,
             CoefSpln *coef)
{
    double f;
    double x;
    int N =size-1;
    double deriv1;
    double deriv2;
    double h1;
    double h0 = arZ[1] - arZ[0];
    double A = 0.;
    double C = 2.*h0;
    double B = h0;
    int i;
    int j;
	//vector<double> Alfa(size);
	//vector<double> Beta(size);
	//vector<double> F(size);
    derivat(arZ, arFz, size, &deriv1, &deriv2);

    coef[0].Alfa = 0.;
    coef[0].Beta = 0.;
    f = 6.*((arFz[1]-arFz[0])/h0-deriv1);

    coef[0].F = f;
    coef[1].Alfa = -B/C;
    coef[1].Beta = f/C;

    for(i=1;i<N;i++ )
    {
        h1 = arZ[i+1] - arZ[i];
        coef[i].F = 6.*((arFz[i+1] - arFz[i])/h1 - (arFz[i] - arFz[i-1])/h0);

        A = h0;
        C = 2.*(h1+h0);
        B = h1;
        x = A*coef[i].Alfa+C;
        coef[i+1].Alfa = -B/x;
        coef[i+1].Beta = (coef[i].F - A*coef[i].Beta)/x;
        h0 = h1;
        //std::cout<<" i = "<< i<< std::endl;
    }
    coef[N].F = 6.*(deriv2 - (arFz[N]-arFz[N-1])/h0);
    A = h0;
    C = 2.*h0;
    B = 0.;
    coef[N].M = (coef[N].F - A*coef[N].Beta)/(A*coef[N].Alfa+C);

    for(j=N-1; j>=0; j--)
    {
        coef[j].M = coef[j+1].Alfa*coef[j+1].M + coef[j+1].Beta;
    }

}

int cspline(const double *arZ,
             const double *arFz,
             const int size,
             CoefSpln *coef)
{
    double h;
    int N = size-1;
    int i;
    double M0;
    if(size<4)
        return 0;
    moment(arZ, arFz, size, coef);

    for(i=0; i<N; i++)
    { 
        h = arZ[i+1]-arZ[i];
        M0 = coef[i].M;
        coef[i].A = (coef[i+1].M - M0)/(6.*h);
        coef[i].B = M0/2.;
		coef[i].C = (arFz[i+1] - arFz[i])/h - coef[i+1].M*h/6. -  M0*h/3.;
        coef[i].D = arFz[i];
        coef[i].Z = arZ[i];
    }
    coef[N].Z = arZ[N];
    return 1;
}


//void citerpForZdef(Motion *motion, CoefSpln *coefX, CoefSpln *coefY)
//{
//    int sizeSpln = motion->sizeI;
//    double xi = motion->XbetaI2[sizeSpln-1];
//    double yi = motion->YbetaI2[sizeSpln-1];
//    int n = 0;
//    double x;
//    int size = motion->size;
//    int i = 0;
//    int j = 0;
//    cspline(motion->ZI,motion->YbetaI2,sizeSpln,coefY);
//    cspline(motion->ZI,motion->XbetaI2,sizeSpln,coefX);
//
//    for( i = 0; i<sizeSpln-1; i++)
//    {
//
//        for( j=0; j<3; j++)
//        {
//            n = i*3 + j + 1 ;
//            x = motion->Z[n] - motion->ZI[i];
//            motion->XbetaI2[n] = coefX->A[i]*x*x*x + coefX->B[i]*x*x
//                                       + coefX->C[i]*x + coefX->D[i];
//
//            motion->YbetaI2[n] = coefY->A[i]*x*x*x + coefY->B[i]*x*x
//                                       + coefY->C[i]*x + coefY->D[i];
//        }
//    }
//
//    motion->XbetaI2[size-1] = xi;
//    motion->YbetaI2[size-1] = yi;
//}
//
