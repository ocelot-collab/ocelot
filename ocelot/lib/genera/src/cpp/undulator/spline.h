#ifndef SPLINE_H
#define SPLINE_H

struct CoefSpln
{
    double A;
    double B;
    double C;
    double D;
    double Z;
    double Alfa; //?
    double Beta; //?
    double F;// ?
    double M;
};

void derivat(const double *arZ,
             const double *arFz,
             const int size,
             double *deriv1,
             double *deriv2);

void moment( const double *arZ,
             const double *arFz,
             const int size,
             CoefSpln coef);

int cspline(const double *arZ,
             const double *arFz,
             const int size,
             CoefSpln *coef);

//void citerpForZdef(Motion *motion, CoefSpln *coefX, CoefSpln *coefY);

//void progon(const TermMatrix &abc, const QVector<double> &F, CoefProgon *coef);

#endif // SPLINE_H
