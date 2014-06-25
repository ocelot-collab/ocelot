#ifndef TRAJECT_H
#define TRAJECT_H

struct ElemMotion
{
    double X;
    double Xbeta;
    double Y;
    double Ybeta;
    double Z;
    double Zbeta;
    double Bx;
    double By;
    double Bz;
};
//#ifndef defined DLL_EXPORT
//#define DECLDIR __declspec(dllexport)
//#else
//#define DECLDIR __declspec(dllimport)
//#endif

int trajectory(double *aMagField,
               const int colsMF,
               const int lenMF,
               const double *misalign,
               const int bRough,
               const double *aInitCond,
               const double *undul_param,
               const int NstepMotion,
               const int Nsuperperiod,
               double *aMotion);


#endif
