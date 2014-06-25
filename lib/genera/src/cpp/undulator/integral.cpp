
#include "integral.h"

int gaussBetaSquare(Motion *motion)
{
    int N = motion->size;
    int size = (N + 1)/3;
    double f1;
    double f2;
    double integ1 = 0.;
    double integ2 = 0.;
    double h2 = (motion->Z[N-1]-motion->Z[0])/2./(size-1);
    double w[] = {0.5555555555555556*h2, 0.8888888888888888*h2, 0.5555555555555556*h2};
    int n;
    int i = 0;
    int j = 0;
    if(motion->bRough>0)
        return 0;
    //mtnParam.sizeI = size;
    motion->XbetaI2[0] = 0.;
    motion->YbetaI2[0] = 0.;
    motion->ZI[0] = motion->Z[0];
    for(i=0; i<size-1;i++)
    {
        for(j=0;j<3; j++)
        {
            n = i*3 + j +1;
            f1 = motion->Xbeta[n];
            f2 = motion->Ybeta[n];
            integ1 += f1*f1*w[j];
            integ2 += f2*f2*w[j];
        }
        motion->XbetaI2[i+1] =  integ1;
        motion->YbetaI2[i+1] =  integ2;
        motion->ZI[i+1] = motion->ZI[i] + h2*2.;
    }
    return 1;
}

void citerpForZdef(Motion *motion)
{
    int sizeSpln = (motion->size + 1)/3;
	motion->ZI = new double[ sizeSpln];
	gaussBetaSquare(motion);

	
    double xi = motion->XbetaI2[sizeSpln-1];
    double yi = motion->YbetaI2[sizeSpln-1];
    int n = 0;
    double x;
    int size = motion->size;
    int i = 0;
    int j = 0;
    
	CoefSpln *coefX = new CoefSpln[sizeSpln];
	CoefSpln *coefY = new CoefSpln[sizeSpln];
    
    cspline(motion->ZI,motion->YbetaI2,sizeSpln,coefY);
    cspline(motion->ZI,motion->XbetaI2,sizeSpln,coefX);

    for( i = 0; i<sizeSpln-1; i++)
    {
        for( j=0; j<3; j++)
        {
            n = i*3 + j + 1 ;
            x = motion->Z[n] - motion->ZI[i];
            motion->XbetaI2[n] = ((coefX[i].A*x + coefX[i].B)*x + coefX[i].C)*x + coefX[i].D;
            motion->YbetaI2[n] = ((coefY[i].A*x + coefY[i].B)*x + coefY[i].C)*x + coefY[i].D;
        }
    }

    motion->XbetaI2[size-1] = xi;
    motion->YbetaI2[size-1] = yi;
	delete []motion->ZI;
	delete []coefX;
	delete []coefY;
}

