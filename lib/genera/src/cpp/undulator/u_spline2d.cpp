#include "u_spline2d.h"
#include <math.h>
//#include "magfield.h"

void searchLongInt1(Field *field, double z)
{
    
	if( z < field->Z[field->nz-2])
		for(int i = field->numIntZ; i<field->nz; i++)
		{
			if (z>= field->Z[i] && z < field->Z[i+1])
			{
				field->numIntZ = i;
				break;
			}
		}
		//while (!(z>= field->Z[i] && z < field->Z[i+1]))
		//{
		//	i++;
		//}
        //#print z, np.argwhere((z - field.Z) >= -1e-12)
        //field.numIntZ = np.argwhere((z - field.Z) >= -1e-12)[-1][0]
    else
        field->numIntZ = field->nz-2;
};

int argwhere(double *arrX, int size, double X)
{
	int i = 0;
	while(!( X >= arrX[i] && X < arrX[i+1]))
		i++;
	return i;
}
void searchTransInt(Field *field, double x, double y, double *k, int *indx)
{
    int bSpline = 0;
    int nx;
	nx = field->numIntX;
    int ny;
	ny = field->numIntY;    
	int ky = 1;
    int kx = 1;  
    if(!(field->X[field->numIntX+1] > x && x >=field->X[field->numIntX]) || 
        !(field->Y[field->numIntY+1] > y && y >=field->Y[field->numIntY]))
	{
        
		ny = argwhere(field->Y, field->ny, y);
        nx = argwhere(field->X, field->nx, x);
        bSpline = 1;    
	}
	
    if (y == field->Y[field->ny - 1])
	{
        ny = ny - 1;
        bSpline = 1;
	}
    if (x == field->X[field->nx - 1])
	{
        nx =  nx - 1;
        bSpline = 1;
	}
    field->numIntX = nx;
    field->numIntY = ny;
  
    double q = (y - field->Y[ny])/(field->Y[ny+1] - field->Y[ny]);
        if((ny <= 0 || q > 0.5) && ny != field->ny - 2)
	{
        ny = ny + 1;
        ky = -1;
        q = 1. - q;
	}
    double p = (x - field->X[nx])/(field->X[nx+1] - field->X[nx]);
        if((nx<=0 || p > 0.5) && nx != field->nx-2)
	{
        nx = nx + 1;
        kx = -1;
        p = 1. - p;
	}
    //# coefficients 
	
    k[0] = q*(q - 1.)/2.;
    k[1] = p*(p - 1.)/2.;
    k[2] = (1. + p*q - p*p - q*q);
    k[3] = p*(p - 2.*q + 1.)/2.;
    k[4] = q*(q - 2.*p + 1.)/2.;
    k[5] = p*q;
    //#indixes 

    indx[0] = nx*field->ny      + ny-ky;
	indx[1] = (nx-kx)*field->ny + ny;
	indx[2] = nx*field->ny      + ny;
	indx[3] = (nx+kx)*field->ny + ny;
	indx[4] = nx*field->ny      + ny+ky;
	indx[5] = (nx+kx)*field->ny + ny+ky;
    //iy = [ny-ky, ny, ny, ny, ny+ky, ny+ky];
    if( ky == -1 || kx == -1 || bSpline == 1)
        splGet(field, indx);
    //return 1;
	
};

void splGet(Field *field, int *indx)
{

    for(int i=0; i < 6; i++) // in range(len(ix)):
	{
		
		if(field->splPack[indx[i]].check == 0) //(  splX[ix[i]*Ny + iy[i]] == 0)
		{
			// allocation memory !!! 
			field->splPack[indx[i]].check = 1;
			field->splPack[indx[i]].splX = new CoefSpln[field->nz];
			field->splPack[indx[i]].splY = new CoefSpln[field->nz];
			field->splPack[indx[i]].splZ = new CoefSpln[field->nz];

			cspline(field->Z, field->magField[indx[i]].Bx, field->nz, field->splPack[indx[i]].splX);
			cspline(field->Z, field->magField[indx[i]].By, field->nz, field->splPack[indx[i]].splY);
			cspline(field->Z, field->magField[indx[i]].Bz, field->nz, field->splPack[indx[i]].splZ);
			
		}
	}
}
    
int spline2d(Field *field, double x, double y, double z, double *Bx, double *By, double *Bz)
{

    //double cos_t = cos(field->undul.tilt);
    //double sin_t = sin(field->undul.tilt);

    if( x>field->X[field->nx-1] || x<field->X[0] || y > field->Y[field->ny-1] || y< field->Y[0])
    {
        //print "x or y out of bounds. spline.clspline2d_py.py"
        return -200;
    }

    double k[6];
    int indx[6];// = new int [6];


    searchTransInt(field, x, y,  k,  indx);// #search indexes 

    if (z > field->Z[field->numIntZ + 1])
        searchLongInt1(field, z);

    int nz = field->numIntZ;
    double h = z - field->Z[nz];

    double h2 = h*h;
    double BX[6];
    double BY[6];
    double BZ[6]; 
    //int Ny = field->ny;
    CoefSpln splx, sply, splz;


    for( int i = 0; i<6; i++) // in range(len(ix)):
    {

        splx = field->splPack[indx[i]].splX[nz];
        sply = field->splPack[indx[i]].splY[nz];
        splz = field->splPack[indx[i]].splZ[nz];

        BX[i] = splx.A*h2*h + splx.B*h2 + splx.C*h + splx.D;
        BY[i] = sply.A*h2*h + sply.B*h2 + sply.C*h + sply.D;
        BZ[i] = splz.A*h2*h + splz.B*h2 + splz.C*h + splz.D;
    }
    /* */
    //#print "spline 2D ", field.X[nx], field.Y[ny]
    *Bx = k[0]*BX[0] + k[1]*BX[1] + k[2]*BX[2] + k[3]*BX[3] + k[4]*BX[4] + k[5]*BX[5];
    *By = k[0]*BY[0] + k[1]*BY[1] + k[2]*BY[2] + k[3]*BY[3] + k[4]*BY[4] + k[5]*BY[5];
    *Bz = k[0]*BZ[0] + k[1]*BZ[1] + k[2]*BZ[2] + k[3]*BZ[3] + k[4]*BZ[4] + k[5]*BZ[5];

    //*Bx = Bx_new*cos_t - By_new*sin_t;
    //*By = Bx_new*sin_t + By_new*cos_t;

    return 1;
}
