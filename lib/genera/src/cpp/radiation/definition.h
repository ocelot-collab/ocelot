#ifndef DEFINITION_H
#define DEFINITION_H

// exemple for pure C

//typedef struct _CoefSpln
//{
//    double A;
//    double M;
//    //int size;
//}CoefSpln;


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



struct Screen
{
	double eStep;
    double eStart;
    int eNstep;
    double xStep;
    double xStart;
    int xNstep;
    double yStep;
    double yStart;
    int yNstep;
    double Distance;
	double Zstart;
    double lperiod;
    int nperiods;
    int status;
    double *ReEx;
    double *ImEx;
    double *ReEy;
    double *ImEy;
	double *Phase;
};





#endif // DEFINITION_H
