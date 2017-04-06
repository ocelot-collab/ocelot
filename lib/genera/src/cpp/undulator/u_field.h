#ifndef FIELD_H
#define FIELD_H
#include "spline.h"

struct UndulParametr
{
    double ampl_Bx;
    double ampl_By;
    double phase;
    int Nperiod;
    double period;
    double ax;
    double start; //for future
    double end; //for future
    double shift_y;
    double shift_x;
    double vangle;
    double hangle;
    double tilt;
};

struct SplPack
{	
	CoefSpln *splX;
	CoefSpln *splY;
	CoefSpln *splZ;
	int check;
};

struct MagField
{
	double *Bx;
	double *By;
	double *Bz;
};

struct Field
{
    UndulParametr undul;
    double *X;
    double *Y;
    double *Z;
    MagField *magField;
    int nx; // number of points in horizoontal direction for field mesh
    int ny; // number of points in vertical direction for field mesh
    int nz; // number of points in longitudinal direction for field mesh
    int cols;    
    SplPack *splPack;
		
    int numIntX;
    int numIntY;
    int numIntZ;

    int Status; // if 1 is planar; if 2 is spiral; if 3 is arbitrary field;
};

int typeOfField(Field *field,  double *aMagField, const int colsMF, const int lenMF);
int fieldAnalytic(Field *field,  double *aMagField, const int colsMF, const int lenMF);
int fieldPlanar(Field *field,  double *aMagField, const int colsMF, const int lenMF);
int fieldSpiral(Field *field,  double *aMagField, const int colsMF, const int lenMF);
int argwhere(double *arrX, double X);
void iterator(double *arrX, int start, int stop, int step, double *newX);

int fieldArbitrary(Field *field, double *aMagField, const int colsMF, const int lenMF);

#endif // FIELD_H
