//#include <windows.h>
#include "spline.h"
#include "u_rungekutta.h"
#include "u_field.h"
#include "u_magfield.h"
#include "motion.h"
#include "u_spline2d.h"
#include "integral.h"
#include <omp.h>
//#include <math.h>
//#include "traject.h"


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
void memoryForLastElem(ElemMotion *elemMotion, Motion motion, int NstepMotion)
{
    elemMotion->X = motion.X[NstepMotion - 2];
    elemMotion->Y = motion.Y[NstepMotion - 2];
    elemMotion->Z = motion.Z[NstepMotion - 2];
    elemMotion->Xbeta = motion.Xbeta[NstepMotion - 2];
    elemMotion->Ybeta = motion.Ybeta[NstepMotion - 2];
    elemMotion->Zbeta = motion.Zbeta[NstepMotion - 2];
    elemMotion->Bx = motion.Bx[NstepMotion - 2];
    elemMotion->By = motion.By[NstepMotion - 2];
    elemMotion->Bz = motion.Bz[NstepMotion - 2];
}

void removeFirstElem(ElemMotion elemMotion, Motion *motion)
{
    motion->X[0] = elemMotion.X;
    motion->Y[0] = elemMotion.Y;
    motion->Z[0] = elemMotion.Z;
    motion->Xbeta[0] = elemMotion.Xbeta;
    motion->Ybeta[0] = elemMotion.Ybeta;
    motion->Zbeta[0] = elemMotion.Zbeta;
    motion->Bx[0] = elemMotion.Bx;
    motion->By[0] = elemMotion.By;
    motion->Bz[0] = elemMotion.Bz;
}

/*
BOOL APIENTRY DllMain(HMODULE, DWORD, LPVOID)
{
    return TRUE;
}
*/
extern "C"
{
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


int da_undulator(double *aMagField,
                 const int colsMF,
                 const int lenMF,
                 const double *misalign,
                 const int bRough,
                 const double *undul_param,
                 const int NstepMotion,
                 const int Nsuperperiod,
                 double gamma,
                 double *c_particle,
                 int np);

}

// this function is main function

int trajectory(double *aMagField,
               const int colsMF,
               const int lenMF,
               const double *misalign,
               const int bRough,
               const double *aInitCond,
               const double *undul_param,
               const int NstepMotion,
               const int Nsuperperiod,
               double *aMotion)
{

    Particle particle;
    Field field;

    // undul_param = [Bx, By, phase, nPeriods,lenPeriod, ax]
    field.undul.ampl_Bx = undul_param[0]; //ampl_Bx;
    field.undul.ampl_By = undul_param[1]; //ampl_By;
    field.undul.phase = undul_param[2];
    field.undul.Nperiod = undul_param[3];
    field.undul.period = undul_param[4];
    field.undul.ax = undul_param[5];

    field.undul.shift_y = misalign[0];    //  Now it works
    field.undul.shift_x = misalign[1];    // Now it works
    field.undul.vangle = misalign[2];    // maybe I will add this future everywhere !!! NOW IT IS NOT USED
    field.undul.hangle = misalign[3];    // maybe I will add this future everywhere !!! NOW IT IS NOT USED
    if(Nsuperperiod > 1)
    {
        field.undul.vangle = 0.;    // maybe I will add this future everywhere !!! NOW IT IS NOT USED
        field.undul.hangle = 0.;    // maybe I will add this future everywhere !!! NOW IT IS NOT USED
    }
    field.undul.tilt = misalign[4];      //  Now it works

    //
    // this function can define type of undulator, depanding on condition it can define four type of undulator
    // 1. planar undulator with uniform transverse field (By(x,y,z) = cubic interpolation in longitudinal direction as function from Z)
    // 2. planar undulator with nonuniform transverse field, you must define parameter a or ksi kx = 2*pi/ksi, and period kz = 2*pi/period
    // for this two cases you can use data from file or define start and end of field map (optional for future)
    // 3. spiral undulator
    // 4. arbitrary magtetic field
    //
    //
    // function return error, if all OK, return 0
    // if it is unable to detirminate type it return '-100'
    //
    int err = typeOfField(&field, aMagField, colsMF,  lenMF);
    if(!err)
        return err;


    particle.X = aInitCond[0]; // initial horizontal condition [mm]
    particle.Y = aInitCond[1]; // initial vertical condition [mm]
    particle.Z =  aInitCond[2]; // initial longitudinal condition [mm]
    particle.betaX = aInitCond[3]; // initial horizontal angel [rad]
    particle.betaY = aInitCond[4]; // initial vertical angel [rad]
    particle.gamma = aInitCond[5]; // gamma is gamma
    particle.mass = 1.;
    particle.charge = -1;
    particle.curr = 1.; // current, don't use here, only for radiation

    int step_size = NstepMotion;
    if(bRough == 0)
    {										/* this condition is very important, if bRough = 0 appear inner cycle in RK */
        step_size = (NstepMotion+1)/3;      /*  for(i=0;i<3;i++) and we must decrease NstepMotion on 3 times */
    }

    int starPosition;
    ElemMotion elemMotion;

    for(int i = 0; i<Nsuperperiod; i++)
    {

        Motion motion;
        //
        // this function just allocates memory for struct motion from array aMotion,
        // like this 	motion->X = &aMotion[0]; motion->Y = &aMotion[NstepMotion];
        //

        starPosition = i*(NstepMotion-(2-bRough));
        initMotion(&motion, bRough, NstepMotion, Nsuperperiod, starPosition, aMotion);

        err = rkgauss(step_size, particle, &field, &motion);

        if( bRough == 0)
        {
            if( i > 0 )
                removeFirstElem(elemMotion, &motion);

            memoryForLastElem(&elemMotion, motion, NstepMotion);

        }


        field.numIntZ = 0;
        field.numIntX = 0;
        field.numIntY = 0;

        particle.X = motion.X[NstepMotion - 1];			// initial horizontal condition [mm]
        particle.Y = motion.Y[NstepMotion - 1];			// initial vertical condition [mm]
        particle.Z = motion.Z[NstepMotion - 1];			// initial longitudinal condition [mm]
        particle.betaX = motion.Xbeta[NstepMotion - 1]; // initial horizontal angel [rad]
        particle.betaY = motion.Ybeta[NstepMotion - 1]; // initial vertical angel [rad]

    }

    Motion motion;

    starPosition = 0;
    initMotion(&motion, bRough, NstepMotion, Nsuperperiod, starPosition, aMotion);

    motion.size = Nsuperperiod*NstepMotion - (Nsuperperiod - 1)*(2-bRough);
    ////gaussBetaSquare(&motion);

    citerpForZdef(&motion);
	

    if(field.Status==3)
    {
        delete []field.X;
        delete []field.Y;
        delete []field.Z;
    }
    if (field.Status != 4)
    {
        for (int i = 0; i<field.nx*field.ny; i++)
        {
            if(field.Status==3)
            {
                delete []field.magField[i].Bx; //np.zeros((lx,ly,lz));
                delete []field.magField[i].By; //np.zeros((lx,ly,lz));
                delete []field.magField[i].Bz;  //np.zeros((lx,ly,lz));
            }

            if(field.splPack[i].check != 0) //(  splX[ix[i]*Ny + iy[i]] == 0)
            {
                if(field.Status==1)
                    delete []field.splPack[i].splY;
                if(field.Status==2)
                {
                    delete []field.splPack[i].splY;
                    delete []field.splPack[i].splX;
                }
                if(field.Status==3)
                {
                    delete []field.splPack[i].splY;
                    delete []field.splPack[i].splX;
                    delete []field.splPack[i].splZ;
                }
            }
        }
	
	delete []field.magField;
	delete []field.splPack;
    }

	return err; //motion.X[15000]*10;
}


int da_undulator(double *aMagField,
                 const int colsMF,
                 const int lenMF,
                 const double *misalign,
                 const int bRough,
                 const double *undul_param,
                 const int NstepMotion,
                 const int Nsuperperiod,
                 double gamma,
                 double *c_particle,
                 int np)
{
    /* Python class Particle
    class Particle:
        def __init__(self, x=0.0, y=0.0, px=0.0, py=0.0, s=0.0, p=0.0,  tau=0.0):
            self.x = x
            self.y = y
            self.px = px       # horizontal (generalized) momentum
            self.py = py       # vertical (generalized) momentum
            self.p = p         # longitudinal momentum
            self.s = s
            self.tau = tau     # time-like coordinate wrt reference particle in the bunch
    */
    int ret = 0;
    #pragma omp parallel for
    for(int i = 0; i < np; i++)
    {
        double x = c_particle[i*7 + 0]; // initial horizontal condition [m]
        double y = c_particle[i*7 + 1]; // initial vertical condition [m]
        double px = c_particle[i*7 + 2]; // initial horizontal angel [rad]
        double py = c_particle[i*7 + 3]; // initial vertical angel [rad]
        double p = c_particle[i*7 + 4]; // momentum
        double s =  c_particle[i*7 + 5]; // initial longitudinal condition [m]
        double tau = c_particle[i*7 + 6];

        double aInitCond[6] = {x*1000,y*1000, s*0., px, py, gamma*(1 + p)};

        int npoints_traj = NstepMotion*Nsuperperiod - (Nsuperperiod - 1)*(2-bRough);
        int nstep = npoints_traj*11;
        double *aMotion = new double [nstep];

        ret = trajectory(aMagField, colsMF, lenMF, misalign, bRough, aInitCond, undul_param, NstepMotion, Nsuperperiod, aMotion);
        //if (ret<0)
        //{
        //   return ret;
        //}

        /* x */   c_particle[i*7 + 0] = aMotion[npoints_traj- 1 ]/1000.; // initial horizontal condition [m]
        /* y */   c_particle[i*7 + 1] = aMotion[2*npoints_traj - 1]/1000.; // initial vertical condition [m]
        /* px */  c_particle[i*7 + 2] = aMotion[4*npoints_traj  - 1]; // initial horizontal angel [rad]
        /* py */  c_particle[i*7 + 3] = aMotion[5*npoints_traj  - 1]; // initial vertical angel [rad]
        /* p */   c_particle[i*7 + 4] = p; // momentum
        /* s */   c_particle[i*7 + 5] = aMotion[3*npoints_traj - 1]/1000.; // initial longitudinal condition [m]
        /* tau */ c_particle[i*7 + 6] = tau;
        delete []aMotion;
    }
    return ret;
}




