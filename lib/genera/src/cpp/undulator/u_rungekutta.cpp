#include "u_rungekutta.h"
#include <math.h>

int rkgauss( const int NstepMotion,
             const Particle particle,
             Field *field,
             Motion *motion)
{
    double X, Y, Z;
    double bx, by, bz;
    double kx1, ky1;
    double kx2, ky2;
    double kx3, ky3;
    double kx4, ky4;
    double mx1, my1;
    double mx2, my2;
    double mx3, my3;
    double mx4, my4;
    double bxconst;
    double byconst;
    int err;
    //double bzconst;
    double dz1;
    double sqrt35 = 0.5*sqrt(3./5.);
    double kz[] = {0.5-sqrt35, sqrt35, sqrt35, 0.5-sqrt35};
    double cmm = 299792.458;
    double k;
    double massElectron = 0.510998910e+6; // rest mass of electron
    double dzk ;
    double Bx;
    double By;
    double Bz;
    double *pBx = &Bx;
    double *pBy = &By;
    double *pBz = &Bz;
    double bx2;
    double by2;
    //int lastNumb=0;
    double dGamma2 = 1. - 0.5/(particle.gamma*particle.gamma);
    double sq;
    int i = -1;
    //double dZperBz;
    double dz;
    int Nrough;
    int n = 0;
    int j = 0;
    int dn = 0;
    k = particle.charge*cmm/(massElectron*particle.mass*particle.gamma);
    double endZ = field->undul.end;
    double startZ = field->undul.start;
    dz1 = (endZ - startZ)/(NstepMotion-1.);

    motion->Xbeta[0] = particle.betaX;// particle.betaX in mrad motion.Xbeta in rad !!!  IT SEEMS than the statement IS WRONG
    motion->Ybeta[0] = particle.betaY;
    motion->Zbeta[0] = dGamma2 - (particle.betaX*particle.betaX + particle.betaY*particle.betaY)/2.;

    motion->X[0] = particle.X;
    motion->Y[0] = particle.Y;
    motion->Z[0] = particle.Z;
    double Zshift = startZ - particle.Z;
	
    if(motion->bRough > 0)
    {
        Nrough = 1;
        kz[0] = 1.;
        dn = 1;
    }
    else
    {
        Nrough = 3;
        dn = 0;
    }
    //motion->XbetaI2 = new double [10000];
    for( n = 0; n < NstepMotion-dn; n++)
    {
        //B[2];
        //std::cout<<"i = "<<i<<std::endl;

        if(n==NstepMotion-1 && motion->bRough<=0)
        {
            kz[0] = kz[0]/2.;
            Nrough = 1;
        }
        for(j=0; j<Nrough; j++)
        {

            dz = dz1*kz[j];
            i++;

            dzk = dz*k;
            X = motion->X[i];
            Y = motion->Y[i];
            Z = motion->Z[i];

            bxconst = motion->Xbeta[i];
            byconst = motion->Ybeta[i];
            bz = motion->Zbeta[i];

            bx = bxconst;
            by = byconst;

            err = B3D(field, X, Y, Z + Zshift, pBx, pBy, pBz);
            if(err<0)
                return err;
            motion->Bx[i] = Bx;
            motion->By[i] = By;
            motion->Bz[i] = Bz;

            //motion->XbetaI2[i] = motion->By[i];
            kx1 = bx*dz;
            ky1 = by*dz;
            bx2 = bx*bx;
            by2 = by*by;
            //sq = 1 + (bx2 + by2)/2.;
            sq = sqrt(1 + bx2 + by2);
            //mx1 = sq*((by*Bz - By)*(1+bx2) + bx*by*(Bx - bx*Bz))*dzk;
            //my1 = sq*((Bx - bx*Bz)*(1+by2) + bx*by*(by*Bz - By))*dzk;
            // the same but shorter exp.
            mx1 = sq*(by*Bz - By*(1+bx2) + bx*by*Bx)*dzk;
            my1 = -sq*(bx*Bz - Bx*(1+by2) + bx*by*By)*dzk;
            //K2

            err = B3D(field, X + kx1/2., Y + ky1/2., Z + dz/2. + Zshift, pBx, pBy, pBz);
            if(err<0)
                return err;
            bx = bxconst + mx1/2.;
            by = byconst + my1/2.;

            kx2 = bx*dz;
            ky2 = by*dz;
            bx2 = bx*bx;
            by2 = by*by;
            //sq = 1 + (bx2 + by2)/2.;
            sq = sqrt(1 + bx*bx + by*by);
            mx2 = sq*(by*Bz - By*(1+bx2) + bx*by*Bx)*dzk;
            my2 = -sq*(bx*Bz - Bx*(1+by2) + bx*by*By)*dzk;
            // K3

            err = B3D(field, X + kx2/2., Y + ky2/2., Z + dz/2. + Zshift, pBx, pBy, pBz);
            if(err<0)
                return err;
            bx = bxconst + mx2/2.;
            by = byconst + my2/2.;

            kx3 = bx*dz;
            ky3 = by*dz;
            bx2 = bx*bx;
            by2 = by*by;
            //sq = 1 + (bx2 + by2)/2.;
            sq = sqrt(1 + bx*bx + by*by);
            mx3 = sq*(by*Bz - By*(1+bx2) + bx*by*Bx)*dzk;
            my3 = -sq*(bx*Bz - Bx*(1+by2) + bx*by*By)*dzk;

            //K4
            err = B3D(field, X + kx3, Y + ky3, Z + dz + Zshift,pBx, pBy, pBz);
            if(err<0)
                return err;
            bx = bxconst + mx3;
            by = byconst + my3;


            kx4 = bx*dz;
            ky4 = by*dz;
            bx2 = bx*bx;
            by2 = by*by;
            //sq = 1 + (bx2 + by2)/2.;
            sq = sqrt(1 + bx*bx + by*by);
            mx4 = sq*(by*Bz - By*(1+bx2) + bx*by*Bx)*dzk;
            my4 = -sq*(bx*Bz - Bx*(1+by2) + bx*by*By)*dzk;
            //            if( n == 0)
            //            {
            //                mexPrintf("SPIRAL UNULATOR, %g.\n", motion->Zbeta[0]);
            //            }
            motion->X[i+1] = X + 1/6.*(kx1 + 2.*kx2 + 2.*kx3 + kx4);
            motion->Xbeta[i+1] = (bxconst + 1/6.*(mx1 + 2.*mx2 + 2.*mx3 + mx4)); // conversion in mrad
            motion->Y[i+1] = Y + 1/6.*(ky1 + 2.*ky2 + 2.*ky3 + ky4);
            motion->Ybeta[i+1] = (byconst + 1/6.*(my1 + 2.*my2 + 2.*my3 + my4));
            motion->Z[i+1] = Z + dz;

            motion->Zbeta[i+1] = dGamma2 - (motion->Xbeta[i+1]*motion->Xbeta[i+1] + motion->Ybeta[i+1]*motion->Ybeta[i+1])/2.; //bz;
            // very strange situation, here was just "  motion->Zbeta[i+1] = bz; " !!!! I think it very big mistake

        }
        if(n==0 && motion->bRough<=0)
        {
            kz[0] = kz[0]*2.;
        }
        //motion->X[i] = dz;
    }
    err = B3D(field, motion->X[i+1], motion->Y[i+1], motion->Z[i+1] + Zshift, pBx, pBy, pBz);
    if(err<0)
        return err;
    //std::cout<<" B = "<<megaPack->motion.Z[i+1] <<std::endl;
    motion->Bx[i+1] = Bx;
    motion->By[i+1] = By;
    motion->Bz[i+1] = Bz;
    //delete []coefY;
    //delete []coefX;
	

    return err;
}

