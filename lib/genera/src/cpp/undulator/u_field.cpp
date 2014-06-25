#include "u_field.h"

int typeOfField(Field *field, double *aMagField, const int colsMF, const int lenMF)
{
    switch (colsMF)
    {
    case 2:
        fieldPlanar(field, aMagField, colsMF, lenMF);
        break;
    case 3:
        fieldSpiral(field, aMagField, colsMF, lenMF);
        break;
    case 6:
        fieldArbitrary(field, aMagField, colsMF, lenMF);
        break;
    case 0:
        fieldAnalytic(field, aMagField, colsMF, lenMF);
        break;
    default:
        return -100;
    }
    return 1;
}

int testSort(double *arrX, int len)
{
    double tmp = arrX[0];
    for(int i = 1; i<len; i++)
    {
        if( arrX[i]>tmp)
            tmp = arrX[i];
        else
            return -300;
    }
    return 1;
}

int fieldAnalytic(Field *field, double *aMagField, const int colsMF, const int lenMF)
{
    field->Status = 4;// # planar undulator
    field->undul.start = 0.;
    field->undul.end = (field->undul.Nperiod + 0.5)*field->undul.period; //remove 0 after checking! but I found that ending poles
    // contribute some error. it we must understand
    if (lenMF == 13)
    {
        field->Status = 13;// # planar undulator
        field->undul.start = 0.;
        field->undul.end = field->undul.period;
    }
    return 1;
}

int fieldPlanar(Field *field, double *aMagField, const int colsMF, const int lenMF)
{
    field->Status = 1;// # planar undulator
    field->Z = &aMagField[0];

    field->magField = new MagField [1];

    field->magField[0].By = &aMagField[lenMF];
    field->numIntZ = 0;
    field->nz = lenMF;
    field->nx = 1;
    field->ny = 1;
    field->cols = colsMF;

    field->splPack = new SplPack [1];
    field->splPack[0].check = 1;
    field->splPack[0].splY = new CoefSpln[field->nz];
    cspline(field->Z, field->magField[0].By,field->nz, field->splPack[0].splY);

    field->undul.start = field->Z[0];
    field->undul.end = field->Z[field->nz-1];

    int err = testSort(field->Z,field->nz);
    return err;
}


int fieldSpiral(Field *field, double *aMagField, const int colsMF, const int lenMF)
{
    field->Status = 2;// # spiral undulator
    field->Z = &aMagField[0];

    field->magField = new MagField [1];
    field->magField[0].By = &aMagField[lenMF];
    field->magField[0].Bx = &aMagField[2*lenMF];
    field->numIntZ = 0;
    field->nz = lenMF;
    field->nx = 1;
    field->ny = 1;
    field->cols = colsMF;

    field->splPack = new SplPack [1];
    field->splPack[0].check = 1;
    field->splPack[0].splY = new CoefSpln[field->nz];

    cspline(field->Z, field->magField[0].By ,field->nz,field->splPack[0].splY);
    field->splPack[0].splX = new CoefSpln[field->nz];
    cspline(field->Z, field->magField[0].Bx ,field->nz,field->splPack[0].splX);

    field->undul.start = field->Z[0];
    field->undul.end = field->Z[field->nz-1];

    int err = testSort(field->Z,field->nz);
    return err;
}

//you must use this function only when you sure that X in [arrX[0], arrX[-1]]
int argwhere(double *arrX, double X)
{
	int i = 0;
	while(!(arrX[i] != X))
		i++;
	return i;
}
void iterator(double *arrX, int start, int stop, int step, double *newX)
{
	int n = 0;
	for(int i = start; i<stop; i = i + step)
	{
		newX[n] = arrX[i];
		n++;
	} 
}
int fieldArbitrary(Field *field,  double *aMagField, const int colsMF, const int lenMF)
{
	//# mb distribute all field data on 2D array of vectors  
	//# for instance Bx[i,j], By[i,j], Bz[i,j],  where i = len(X) and j = len(Y) 
	//# and By[i,j] - vector length of len(Z) 
	//# there are three variants of placement x y z
	field->Status = 3;// # arbitrary field 

	field->numIntX = 0;
	field->numIntY = 0;
	field->numIntZ = 0;
	int N = lenMF;
	double *xx = &aMagField[0];
	double *yy = &aMagField[lenMF];
	double *zz = &aMagField[2*lenMF];
	int ki = argwhere(xx, xx[0]); //argwhere(xx!=xx[0]);
	int kj  = argwhere(yy, yy[0]);
	int kk = argwhere(zz, zz[0]);
	//#print ki, kj, kk
	
	if (ki < kj && kj < kk)
	{
		field->nx = kj;
        field->ny = kk/kj;
		field->nz = N/kk;
		field->X = new double [field->nx];
		field->Y = new double [field->ny];
		field->Z = new double [field->nz];
		iterator( xx, 0, kj, 1, field->X); //xx[0:kj]; //#x[0:lx]
		iterator( yy, 0, kk, kj, field->Y); //[0:kk:kj]
		iterator( zz, 0, N, kk, field->Z); //[0:N:kk]
	}
	if (kk < kj && kj < ki)
	{
		field->nx = N/kk;
		field->ny = kk/kj;
		field->nz = kj;
		field->X = new double [field->nx];
		field->Y = new double [field->ny];
		field->Z = new double [field->nz];
		iterator( xx, 0, N, kk, field->X); //xx[0:N:kk];
		iterator( yy, 0, kk, kj, field->Y); //yy[0:kk:kj];
		iterator( zz, 0, kj, 1, field->Z); //zz[0:kj] #z[0:lz];

	}

	int ntrans = field->nx*field->ny; 
	field->splPack = new SplPack [ntrans];

	field->magField = new MagField[ntrans];

	for (int i = 0; i<ntrans; i++)
	{
		field->magField[i].Bx = new double [field->nz]; //np.zeros((lx,ly,lz));
		field->magField[i].By = new double [field->nz]; //np.zeros((lx,ly,lz));
		field->magField[i].Bz = new double [field->nz]; //np.zeros((lx,ly,lz));
		field->splPack[i].check = 0;
	}
	
	for (int i = 0; i < field->nx; i++) //( i in range(lx))
		for(int j = 0; j < field->ny; j++) // j in range(ly)
		{
			int index = ki*i + kj*j;
			int n = 0;
			for (int k = index; k<N; k = k+kk)
			{
				 field->magField[i*field->ny+j].Bx[n] = aMagField[3*lenMF + k];
				 field->magField[i*field->ny+j].By[n] = aMagField[4*lenMF + k];
				 field->magField[i*field->ny+j].Bz[n] = aMagField[5*lenMF + k];
				 n++;
			}
		}


    field->undul.start = field->Z[0];
	field->undul.end = field->Z[field->nz-1];


	int err = testSort(field->Z,field->nz);
	if(err < 0)
		return err;
	err = testSort(field->X,field->nx);
	if(err < 0)
		return err-1;
	err = testSort(field->Y,field->ny);
	if(err < 0)
		return err-2;
	return 1;
}
