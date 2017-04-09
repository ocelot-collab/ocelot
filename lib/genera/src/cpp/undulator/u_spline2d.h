#ifndef SPLINE2D_H
#define SPLINE2D_H

#include "spline.h"
#include "u_field.h"
#include "u_magfield.h"

//struct Field;
int argwhere(double *arrX, int size, double X);
void searchTransInt(Field *field, double x, double y, double *k, int *indx);
void splGet(Field *field, int *indx);
int spline2d(Field *field, double x, double y, double z, double *Bx, double *By, double *Bz);
void searchLongInt1(Field *field, double z);
#endif // SPLINE2D_H