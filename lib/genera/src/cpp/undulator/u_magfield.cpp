#include "u_magfield.h"
#include <math.h>
//#include "iostream"

void searchLongInt(Field *field,  double Z)
{
    int i;
    i = field->numIntZ;
    int size = field->nz;
    if (Z < field->Z[size-2])
        while ((Z > field->Z[i]))
        {
            i = i+1;
            field->numIntZ = i;
        }
    else
        field->numIntZ = size-2;
}


int mixB3D(Field *field, const double x, const double y, const double z,
             double *Bx, double *By, double *Bz)
{
    double dz;
    double Bzin;
    double Bzdin;
    double Bxin;
    double pi = 3.14159265358979;
    double Lu = field->undul.period;
    double By0 = field->undul.ampl_By;
    double Bx0 = field->undul.ampl_Bx;
    //double ksi = field->undul.ax;
    double kz = 2.*pi/Lu;
    double kx = 1/field->undul.ax;
    double ky = sqrt(kz*kz + kx*kx);
    double cosx = cos(kx*x);
    //double ky2 = kz*kz + kx*kx;
    double sinhy = sinh(ky*y);
    double dz2;
    //double sinz = sin(kz*z);
    int i;
    CoefSpln spl;


    //double cos_t = cos(field->undul.tilt);
    //double sin_t = sin(field->undul.tilt);


    if(field->Status == 2 || field->Status == 1) // spiral or planar undulator
    {

        if (z > field->Z[field->numIntZ + 1])
            searchLongInt(field, z);

        i = field->numIntZ;
        spl = field->splPack[0].splY[i];
        dz = z - spl.Z;
        dz2 = dz*dz;
        Bzin = spl.A*dz2*dz + spl.B*dz2 + spl.C*dz + spl.D;
        Bzdin = 3.*spl.A*dz*dz + 2.*spl.B*dz + spl.C;
        *Bx = - kx/ky*sin(kx*x)*sinhy*Bzin;                //
        *By = cosx*cosh(ky*y)*Bzin;          // approximate expression
        *Bz = cosx*sinhy*Bzdin/ky;

        //double Bx_new = *Bx*cos_t - *By*sin_t;
        //double By_new = *Bx*sin_t + *By*cos_t;
        //*Bx = Bx_new*cos_t - By_new*sin_t;
        //*By = Bx_new*sin_t + By_new*cos_t;

        CoefSpln splX;
        if(field->Status == 2)
        {
            splX = field->splPack[0].splX[i];
            Bxin = splX.A*dz2*dz + splX.B*dz2 + splX.C*dz + splX.D;
            *Bx = Bxin;
            *By = Bzin;					// approximate expression
            //*Bx = Bxin*cos_t - Bzin*sin_t;
            //*By = Bxin*sin_t + Bzin*cos_t;
            *Bz = 0.;
            return 1;
        }
        if(field->undul.ax <= 0)
        {
            *Bx = 0.;//- Bzin*sin_t;
            *By = Bzin; //*cos_t;;
            *Bz = 0.;

        }
        return 1;
    }
    if (field->Status == 4 && field->undul.ax > 0) //analytic field
    {
        double sinhy = sinh(ky*y);
        //double cosz = cos(kz*z);
        double strt = 0;
        double end = (field->undul.Nperiod + 0.5)*Lu;
        // delete this after checking and remove 0 in fieldAnalytic() and uncomment below
        //double fld = By0*cos(kz*(z));
        //double fld_1 = -By0*kz*sin(kz*(z));

        double fld = -By0*sin(kz*(z-Lu/2.));
        double fld_1 = -By0*kz*cos(kz*(z-Lu/2.));

        if (z < strt || z > end)
        {
            fld = 0.;
            fld_1 = 0.;
        }
        if (z >= strt && z < strt + Lu/2.)
        {
            fld = By0/2.*sin(kz*z);
            fld_1 = By0/2.*kz*cos(kz*z);
        }
        if (z <= end && z > end - Lu/2.)
        {
            fld = By0/2.*sin(kz*(z - end - Lu/2.));
            fld_1 = By0/2.*kz*cos(kz*(z - end - Lu/2.));
        }

        *Bx = -kx/ky*sin(kx*x)*sinhy*fld;
        *By = cosx*cosh(ky*y)*fld;
        *Bz = 1./ky*cosx*sinhy*fld_1;
        //*Bx = Bx_new*cos_t - By_new*sin_t;
        //*By = Bx_new*sin_t + By_new*cos_t;


        //double cos_alpha = cos(field->undul.vangle);
        //double sin_alpha = sin(field->undul.vangle);
        //*Bx = Bx_t;  //
        //*By = By_t*cos_alpha - Bz_t*sin_alpha;          // exact expression
        //*Bz = Bz_t*cos_alpha + By_t*sin_alpha;       //
    }
    else if (field->Status == 13) // fast calculation for ideal undulator
    {
        double phase = field->undul.phase;
        double fld_y = By0*sin(kz*z);
        double fld_x = Bx0*sin(kz*z + phase);
        //*Bx = fld_x*cos_t - fld_y*sin_t;
        //*By = fld_x*sin_t + fld_y*cos_t;
        *Bx = fld_x;  //
        *By = fld_y;          // exact expression
        *Bz = 0.;       //
    }
    else
    {

        double strt = 0;
        double end = (field->undul.Nperiod + 0.5)*Lu;
        double phase = field->undul.phase;
        // delete this after checking and remove 0 in fieldAnalytic() and uncomment below
        // double fld_y = By0*cos(kz*(z));
        // double fld_x = Bx0*cos(kz*(z) + phase);
        // try to check. uncomment after checking

        double fld_y = -By0*sin(kz*(z-Lu/2.));
        double fld_x = -Bx0*sin(kz*(z-Lu/2.) + phase);
        if (z < strt || z > end)
        {
            fld_y = 0.;
            fld_x = 0;
        }
        if (z >= strt && z < strt + Lu/2.)
        {
            fld_y = By0/2.*sin(kz*z);
            fld_x = Bx0/2.*sin(kz*z + phase);
        }
        if (z <= end && z > end - Lu/2.)
        {
            fld_y = By0/2.*sin(kz*(z - end - Lu/2.));
            fld_x = Bx0/2.*sin(kz*(z - end - Lu/2.) + phase);
        }
        //*Bx = fld_x*cos_t - fld_y*sin_t;
        //*By = fld_x*sin_t + fld_y*cos_t;
        *Bx = fld_x;  //
        *By = fld_y;          // exact expression
        *Bz = 0.;       //
    }
    return 1;

}

int B3D(Field *field, const double x, const double y, const double z,
             double *Bx, double *By, double *Bz)
{
    /*
    Rotation transform for undulator
    x' = x*cos(theta) + y*sin(theta)
    y' = -x*sin(theta) + y*cos(theta)
    */
    double length = field->undul.end - field->undul.start;

    double cos_t = cos(field->undul.tilt);
    double sin_t = sin(field->undul.tilt);

    // title + shift
    double x_t = x*cos_t + y*sin_t - field->undul.shift_x;
    double y_t = -x*sin_t + y*cos_t - field->undul.shift_y;// - field->undul.vangle*z;

    double z_poz = (z - field->undul.start) - length/2.;
    // vertical angle
    double z_v = z + y_t*field->undul.vangle;
    double y_t_v = -z_poz*field->undul.vangle + y_t;
    // horizontal angle
    double z_v_h = z_v + x_t*field->undul.hangle;
    double x_t_h = -z_poz*field->undul.hangle + x_t;
    double Bx1, By1, Bz1;
    int err;
    if (field->Status == 3)
    {
        err = spline2d(field, x_t_h, y_t_v, z_v_h, &Bx1, &By1, &Bz1);
    }
    else
    {
        err = mixB3D(field, x_t_h, y_t_v, z_v_h, &Bx1, &By1, &Bz1);
    }
    // tilt for field
    // rotating vector look like clockwise
    double Bx_t = Bx1*cos_t - By1*sin_t;
    double By_t = Bx1*sin_t + By1*cos_t;
    // vertical
    *By = By_t + Bz1*field->undul.vangle;
    double Bz_v = Bz1 - By_t*field->undul.vangle;
    // horizontal
    *Bx = Bx_t + Bz_v*field->undul.hangle;
    *Bz = Bz_v - Bx_t*field->undul.hangle;
    return err;
}
