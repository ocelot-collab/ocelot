#ifndef MAGFIELD_H
#define MAGFIELD_H
//#include "definition.h"
#include "spline.h"
#include "u_spline2d.h"
#include "u_field.h"



void searchLongInt(Field *field,  double Z);
int mixB3D(Field *field, const double x, const double y, const double z,
             double *Bx, double *By, double *Bz);
int B3D(Field *field, const double x, const double y, const double z,
             double *Bx, double *By, double *Bz);

#endif // MAGFIELD_H
