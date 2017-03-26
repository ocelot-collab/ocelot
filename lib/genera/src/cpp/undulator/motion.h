#ifndef MOTION_H
#define MOTION_H

struct Particle
{
    double X;      /* initial     */
    double Y;
    double Z;      /* coordinates */

    double betaX; // transvers velosity component bettaX = vx/c
    double betaY; // transvers velosity component bettaX = vy/c
    double gamma;
    double mass;   // rest mass in units of electron rest mass: =1 for electron
    int charge;    // charge -1 or +1
    double curr;  // current in mA
};
struct Motion
{
    double *X;
    double *Xbeta;
    double *XbetaI2;
    double *Y;
    double *Ybeta;
    double *YbetaI2;
    double *Z;
    double *Zbeta;
    double *ZI;
    double *Bx;
    double *By;
    double *Bz;
    int size;
    int sizeI;
    int bRough;
};

void initMotion(Motion *motion, const int bRough, const int NstepMotion, int Nsuper, int startPosition, double *aMotion);
#endif // MOTION_H