__author__ = 'Sergey Tomin'

#import numpy as np
from numpy import sqrt, matrix, cos, sin, log, tan, eye, zeros, pi, array, linspace, dot, abs, random, arctan, sign
from scipy import weave
#import os
#os.environ["CC"] = "gcc-4.8"
#os.environ["CXX"] = "gcc-4.8"
def track_und_sym(y0, z, kz, kx, Kx, energy):
    gamma = energy*1957.
    #if energy == 0:
    #    h0 = 0.
    #else:
    h0 = 1./(gamma/Kx/kz)

    #c = 299792458
    #m0 = 0.510998928*1e+6
    #B0 = Kx*m0*kz/c
    #y = array(y0)
    ky = sqrt(kz*kz + kx*kx)

    h = z[1]-z[0]
    N = len(z)
    #print kx, ky, kz, N, 1./h0, h
    #Ax =  1/kz * np.cos(kx*y0[0])*np.cosh(ky*y0[2])*np.sin(kz*z[0])
    #Ay =  kx/(ky*kz) * np.sin(kx*y0[0])*np.sinh(ky*y0[2])*np.sin(kz*z[0])
    q = array([y0[0],y0[1] ,y0[2],y0[3], y0[4], y0[5], kx, ky, kz, h, N, h0])
    #print N
    u = zeros(6)

    code = """
    double x = Q1(0);
    double px = Q1(1);
    double y = Q1(2);
    double py = Q1(3);
    double tau = Q1(4);
    double delta = Q1(5);
    U1(5) = delta;

    double kx = Q1(6);
    double ky = Q1(7);
    double kz = Q1(8);

    double h = Q1(9);
    int N = Q1(10);
    double h0 = Q1(11); //rho = Q1(11);
    double h02 = h0*h0;
    int i;
    double kx2 = kx*kx;
    double ky2 = ky*ky;
    double kz2 = kz*kz;
    //double Kx = kx2/kz2*h02;
    //double Ky = ky2/kz2*h02;
    //double x2, y2, x4, y4;
    h = h/(1. + delta);
    for ( i = 0; i<N-1; i++){
        //x2 = x*x;
        //y2 = y*y;
        //x4 = x2*x2;
        //y4 = y2*y2;
        double chx = cosh(kx*x);
        double chy = cosh(ky*y);
        double shx = sinh(kx*x);
        double shy = sinh(ky*y);
        px = px - h/2.*chx*shx*(kx*ky2*chy*chy + kx2*kx*shy*shy)/(ky2*kz2)*h02;
        py = py - h/2.*chy*shy*(ky2*chx*chx + kx2*shx*shx)/(ky*kz2)*h02;
        //px = px - h*Kx*(x/2. + kx2*x2*x/3. + kz2*x*y2 + (kx2*kx2* x4*x)/15. + (kx2* kz2* x2*x *y2)/3. + (kz2*ky2*x*y4)/6. );
        //py = py - h*Ky*(y/2. + ky2*y2*y/3. + kz2*kx2*x2*y/ky2 + (kx2 *kz2* x2* y2*y)/3. + (kx2*kx2* kz2* x4*y)/(6.*ky2)  + (ky2*ky2* y4*y)/15. );
        //tau = tau - h/(1.+delta) *((px*px + py*py)/2. + Kx*x2/4. + Ky*y2/4.
        //                        + Kx*kx2*x4/12. + Ky*ky2*y4/12.  + Kx*kx2*kx2/90.*x4*x2 + Ky*ky2*ky2/90.*y4*y2
        //                        + h02*(kx2*ky2/12.*x2*y4 + kx2*kx2/12.*x4*y2 + kx2*x2*y2/4.) );
        tau = tau - h/2./(1.+delta) * ((px*px + py*py) + chx*chx*chy*chy/(2.*kz2)*h02 + shx*shx*shy*shy*kx2/(2.*ky2*kz2)*h02);
        x = x + h*px;
        y = y + h*py;

    }
    U1(0) = x;
    U1(1) = px;
    U1(2) = y;
    U1(3) = py;
    U1(4) = tau;
    """
    weave.inline(code, ['u',"q"])
    x = u[0]
    px = u[1]
    y = u[2]
    py = u[3]
    tau = u[4]
    return x, px,y, py, tau

def track_und_mag(y0, z, kz, kx, Kx, energy):
    gamma = energy*1957.
    rho = gamma/Kx/kz
    #c = 299792458
    #m0 = 0.510998928*1e+6
    #B0 = Kx*m0*kz/c

    #y = array(y0)
    ky = sqrt(kz*kz + kx*kx)

    h = z[1]-z[0]
    N = len(z)
    #Ax =  1/kz * np.cos(kx*y0[0])*np.cosh(ky*y0[2])*np.sin(kz*z[0])
    #Ay =  kx/(ky*kz) * np.sin(kx*y0[0])*np.sinh(ky*y0[2])*np.sin(kz*z[0])
    q = array([y0[0],y0[1] ,y0[2],y0[3], y0[4], y0[5], kx, ky, kz, h, N, rho])
    #print N
    u = zeros(6)

    code = """
    double x = Q1(0);
    double px = Q1(1);
    double y = Q1(2);
    double py = Q1(3);
    double tau = Q1(4);
    double delta = Q1(5);
    U1(5) = delta;

    double kx = Q1(6);
    double ky = Q1(7);
    double kz = Q1(8);

    double h = Q1(9);
    int N = Q1(10);
    double rho = Q1(11);

    int i;
    double kx2 = kx*kx;
    double ky2 = ky*ky;
    double kz2 = kz*kz;
    double K = h/(2.*rho*rho*kz2);
    for ( i = 0; i<N-1; i++){

        x = x + h*px;
        y = y + h*py;
        px = px - K* (kx2*x - ky2*ky2/2.*x*y*y + (1./(8.*rho*rho) - kx2/2.)*kz2*x*x*x);
        //py = py - K*(y + kz2*y*y*y/6.);
        py = py - K* (ky2*y + ky2*ky2/6.*y*y*y - (3./(8.*rho*rho) - 3.*kx2/2.)*kz2*x*x*y);
    }
    U1(0) = x;
    U1(1) = px;
    U1(2) = y;
    U1(3) = py;
    U1(4) = tau;
    """
    weave.inline(code, ['u',"q"])
    x = u[0]
    px = u[1]
    y = u[2]
    py = u[3]
    tau = u[4]
    return x, px,y, py, tau


def track_und_RK(y0, z, kz, kx ,Kx, energy):
    gamma = energy*1957.
    #rho = gamma/Kx/kz
    c = 299792458
    m0 = 0.510998928*1e+6
    B0 = Kx*m0*kz/c
    #y = array(y0)
    ky = sqrt(kz*kz - kx*kx)
    #print "B0 = ", B0
    #print "rho = ", rho
    #print "kz = ", kz
    #print "kx = ", kx
    #print "gamma = ", gamma
    h = z[1]-z[0]
    N = len(z)
    #Ax =  1/kz * np.cos(kx*y0[0])*np.cosh(ky*y0[2])*np.sin(kz*z[0])
    #Ay =  kx/(ky*kz) * np.sin(kx*y0[0])*np.sinh(ky*y0[2])*np.sin(kz*z[0])
    q = array([y0[0],y0[1] ,y0[2],y0[3], y0[4], y0[5], kx, ky, kz, h, N, B0, gamma*(1+y0[5])])
    #print N

    u = zeros(N*6)
    code = """
    double charge = 1;
    double mass = 1; //in electron mass

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

    double cmm = 299792458;
    double k;
    double massElectron = 0.510998910e+6; // rest mass of electron
    double dzk ;
    double Bx;
    double By;
    double Bz;

    double bx2;
    double by2;

    double x = Q1(0);
    double px = Q1(1);
    double y = Q1(2);
    double py = Q1(3);
    //double tau = Q1(4);
    //double delta = Q1(5);
    double kx = Q1(6);
    double ky = Q1(7);
    double kz = Q1(8);

    double dz = Q1(9);
    int N = Q1(10);
    double B0 = Q1(11);
    double gamma = Q1(12);
    double z = 0.;
    double dGamma2 = 1. - 0.5/(gamma*gamma);
    double pz = dGamma2 - (px*px + py*py)/2.;

    double sq;
    int i;

    k = charge*cmm/(massElectron*mass*gamma);

    U1(0) = Q1(0);
    U1(1) = Q1(1);
    U1(2) = Q1(2);
    U1(3) = Q1(3);
    U1(4) = z;
    U1(5) = pz;
    dzk = dz*k;
    for(i = 0; i < N-1; i++)
    {


        X = U1(i*6 + 0);
        Y = U1(i*6 + 2);
        Z = U1(i*6 + 4);

        bxconst = U1(i*6 + 1);
        byconst = U1(i*6 + 3);
        bz = U1(i*6 + 5);

        bx = bxconst;
        by = byconst;

        z = Z;
        x = X;
        y = Y;
        Bx = -B0*kx/ky*sin(kx*x)*sinh(ky*y)*cos(kz*z); // here kx is only real
        By = B0*cos(kx*x)*cosh(ky*y)*cos(kz*z);
        Bz = -B0*kz/ky*cos(kx*x)*sinh(ky*y)*sin(kz*z);

        //motion->XbetaI2[i] = motion->By[i];
        kx1 = bx*dz;
        ky1 = by*dz;
        bx2 = bx*bx;
        by2 = by*by;
        //sq = 1 + (bx2 + by2)/2.;
        sq = sqrt(1 + bx2 + by2);

        mx1 = sq*(by*Bz - By*(1+bx2) + bx*by*Bx)*dzk;
        my1 = -sq*(bx*Bz - Bx*(1+by2) + bx*by*By)*dzk;
        //K2
        //err = B3D(field, X + kx1/2., Y + ky1/2., Z + dz/2. + Zshift, pBx, pBy, pBz);
        x = X + kx1/2.;
        y = Y + ky1/2.;
        z = Z + dz/2.;
        Bx = -B0*kx/ky*sin(kx*x)*sinh(ky*y)*cos(kz*z); // here kx is only real
        By = B0*cos(kx*x)*cosh(ky*y)*cos(kz*z);
        Bz = -B0*kz/ky*cos(kx*x)*sinh(ky*y)*sin(kz*z);

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
        //err = B3D(field, X + kx2/2., Y + ky2/2., Z + dz/2. + Zshift, pBx, pBy, pBz);

        x = X + kx2/2.;
        y = Y + ky2/2.;
        z = Z + dz/2.;
        Bx = -B0*kx/ky*sin(kx*x)*sinh(ky*y)*cos(kz*z); // here kx is only real
        By = B0*cos(kx*x)*cosh(ky*y)*cos(kz*z);
        Bz = -B0*kz/ky*cos(kx*x)*sinh(ky*y)*sin(kz*z);

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
        //err = B3D(field, X + kx3, Y + ky3, Z + dz + Zshift,pBx, pBy, pBz);

        x = X + kx3;
        y = Y + ky3;
        z = Z + dz;
        Bx = -B0*kx/ky*sin(kx*x)*sinh(ky*y)*cos(kz*z); // here kx is only real
        By = B0*cos(kx*x)*cosh(ky*y)*cos(kz*z);
        Bz = -B0*kz/ky*cos(kx*x)*sinh(ky*y)*sin(kz*z);

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

        U1((i+1)*6 + 0) = X + 1/6.*(kx1 + 2.*kx2 + 2.*kx3 + kx4);
        U1((i+1)*6 + 1) = bxconst + 1/6.*(mx1 + 2.*mx2 + 2.*mx3 + mx4); // conversion in mrad
        U1((i+1)*6 + 2) = Y + 1/6.*(ky1 + 2.*ky2 + 2.*ky3 + ky4);
        U1((i+1)*6 + 3) = byconst + 1/6.*(my1 + 2.*my2 + 2.*my3 + my4);

        U1((i+1)*6 + 4) = Z + dz;
        U1((i+1)*6 + 5) = dGamma2 - (U1((i+1)*6 + 1)*U1((i+1)*6 + 1) + U1((i+1)*6 + 3)*U1((i+1)*6 + 3))/2.; //bz;


    }
    //err = B3D(field, motion->X[i+1], motion->Y[i+1], motion->Z[i+1] + Zshift, pBx, pBy, pBz);
    //std::cout<<" B = "<<megaPack->motion.Z[i+1] <<std::endl;
    //motion->Bx[i+1] = Bx;
    //motion->By[i+1] = By;
    //motion->Bz[i+1] = Bz;


    """
    weave.inline(code, ['u',"q"])
    x = u[::6]
    y = u[2::6]
    #print x
    #print y
    px = u[1::6]
    py = u[3::6]
    z = u[4::6]
    pz = u[5::6]
    return x, px,y, py, z, pz


def track_und_openmp(u, l, N, kz, kx ,Kx, energy):
    gamma = energy*1957.
    #rho = gamma/Kx/kz
    c = 299792458
    m0 = 0.510998928*1e+6
    B0 = Kx*m0*kz/c

    ky = sqrt(kz*kz - kx*kx)
    z = linspace(0, l, num=N)
    h = z[1]-z[0]
    N = len(z)

    nparticles = len(u)/6
    q = array([ kx, ky, kz, h, N, B0, gamma, nparticles])
    code = """
    double charge = 1;
    double mass = 1; //in electron mass

    double cmm = 299792458;

    double massElectron = 0.510998910e+6; // rest mass of electron

    double kx = Q1(0);
    double ky = Q1(1);
    double kz = Q1(2);

    double dz = Q1(3);
    int N = Q1(4);
    double B0 = Q1(5);

    int npart = Q1(7);

    double sq;
    #pragma omp parallel for
    for(int n = 0; n < npart; n++)
    {
        double kx1, ky1;
        double kx2, ky2;
        double kx3, ky3;
        double kx4, ky4;
        double mx1, my1;
        double mx2, my2;
        double mx3, my3;
        double mx4, my4;
        double bx2;
        double by2;
        double gamma = Q1(6);//*(1+U1(n*6 + 5));
        double dGamma2 = 1. - 0.5/(gamma*gamma);
        double k = charge*cmm/(massElectron*mass*gamma);

        double X = U1(n*6 + 0);
        double bxconst = U1(n*6 + 1);
        double Y = U1(n*6 + 2);
        double byconst = U1(n*6 + 3);
        double Bx;
        double By;
        double Bz;
        //double pz = dGamma2 - (bxconst*bxconst + byconst*byconst)/2.;
        //bz = pz;
        double dzk = dz*k;
        double Z = 0.;

        for(int i = 0; i < N-1; i++)
        {
            double x = X;
            double y = Y;
            double z = Z;
            double bx = bxconst;
            double by = byconst;

            Bx = -B0*kx/ky*sin(kx*x)*sinh(ky*y)*cos(kz*z); // here kx is only real
            By = B0*cos(kx*x)*cosh(ky*y)*cos(kz*z);
            Bz = -B0*kz/ky*cos(kx*x)*sinh(ky*y)*sin(kz*z);

            //motion->XbetaI2[i] = motion->By[i];
            kx1 = bx*dz;
            ky1 = by*dz;
            bx2 = bx*bx;
            by2 = by*by;
            //sq = 1 + (bx2 + by2)/2.;
            sq = sqrt(1 + bx2 + by2);

            mx1 = sq*(by*Bz - By*(1+bx2) + bx*by*Bx)*dzk;
            my1 = -sq*(bx*Bz - Bx*(1+by2) + bx*by*By)*dzk;
            //K2
            //err = B3D(field, X + kx1/2., Y + ky1/2., Z + dz/2. + Zshift, pBx, pBy, pBz);
            x = X + kx1/2.;
            y = Y + ky1/2.;
            z = Z + dz/2.;
            Bx = -B0*kx/ky*sin(kx*x)*sinh(ky*y)*cos(kz*z); // here kx is only real
            By = B0*cos(kx*x)*cosh(ky*y)*cos(kz*z);
            Bz = -B0*kz/ky*cos(kx*x)*sinh(ky*y)*sin(kz*z);

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
            //err = B3D(field, X + kx2/2., Y + ky2/2., Z + dz/2. + Zshift, pBx, pBy, pBz);

            x = X + kx2/2.;
            y = Y + ky2/2.;
            z = Z + dz/2.;
            Bx = -B0*kx/ky*sin(kx*x)*sinh(ky*y)*cos(kz*z); // here kx is only real
            By = B0*cos(kx*x)*cosh(ky*y)*cos(kz*z);
            Bz = -B0*kz/ky*cos(kx*x)*sinh(ky*y)*sin(kz*z);

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
            //err = B3D(field, X + kx3, Y + ky3, Z + dz + Zshift,pBx, pBy, pBz);

            x = X + kx3;
            y = Y + ky3;
            z = Z + dz;
            Bx = -B0*kx/ky*sin(kx*x)*sinh(ky*y)*cos(kz*z); // here kx is only real
            By = B0*cos(kx*x)*cosh(ky*y)*cos(kz*z);
            Bz = -B0*kz/ky*cos(kx*x)*sinh(ky*y)*sin(kz*z);

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

            X = X + 1/6.*(kx1 + 2.*kx2 + 2.*kx3 + kx4);
            bxconst = bxconst + 1/6.*(mx1 + 2.*mx2 + 2.*mx3 + mx4); // conversion in mrad
            Y= Y + 1/6.*(ky1 + 2.*ky2 + 2.*ky3 + ky4);
            byconst = byconst + 1/6.*(my1 + 2.*my2 + 2.*my3 + my4);

            Z = Z + dz;
            //bz = dGamma2 - (bxconst*bxconst +byconst*byconst)/2.; //bz;
        }
        U1(n*6 + 0) = X;
        U1(n*6 + 1) = bxconst;
        U1(n*6 + 2) = Y;
        U1(n*6 + 3) = byconst;
    }



    """
    weave.inline(code, ['u',"q"],
    extra_compile_args =['-O3 -fopenmp'],
    compiler = 'gcc',
    libraries=['gomp'],
    headers=['<omp.h>']
    )
    return u



def track_in_undul_test(y0, z, kz, kx, rho):

    #y = array(y0)
    ky = sqrt(kz*kz + kx*kx)

    h = z[1]-z[0]
    N = len(z)
    Ax =  1/kz * np.cos(kx*y0[0])*np.cosh(ky*y0[2])*np.sin(kz*z[0])
    Ay =  kx/(ky*kz) * np.sin(kx*y0[0])*np.sinh(ky*y0[2])*np.sin(kz*z[0])
    q = np.array([y0[0],y0[1] + Ax/rho ,y0[2],y0[3] + Ay/rho, kx, ky, kz, h, N, rho])
    #print N
    u = zeros(N*4)

    code = """

    U1(0) = Q1(0);
    U1(1) = Q1(1);
    U1(2) = Q1(2);
    U1(3) = Q1(3);

    double kx = Q1(4);
    double ky = Q1(5);
    double kz = Q1(6);

    double h = Q1(7);
    int N = Q1(8);
    double rho = Q1(9);
    int i;
    double z = 0.;
    double h_2 = h/2.;
    double  Ax, Ay, A1, A2;
    double x, y, x1, y1, Px1, Py1, px, py;
    double sinh_y, cos_x, sin_z, cosh_y, sin_x;
    for ( i = 0; i<N-1; i++){

        x = U1(4*i + 0);
        y = U1(4*i + 2);
        px = U1(4*i + 1);
        py = U1(4*i + 3);

        sin_z = sin(kz*z)/rho;
        cosh_y = cosh(ky*y);
        sin_x = sin(kx*x);
        sinh_y = sinh(ky*y);
        cos_x = cos(kx*x);


        Ax = 1./kz*cos_x*cosh_y*sin_z;
        Ay = kx/(ky*kz)*sin_x*sinh_y*sin_z;

        x1 = px - Ax;
        y1 = py - Ay;

        A1 = kx/kz*sin_x*cosh_y*sin_z;
        A2 = cos_x*sinh_y*sin_z;


        Px1 = -x1*A1 + y1*kx*kx/(ky*kz)*A2;
        Py1 = y1*A1 + x1*ky/kz*A2;

        px = U1(4*i + 1) + Px1*h_2;
        py = U1(4*i + 3) + Py1*h_2;
        x = U1(4*i + 0) + x1*h_2;
        y = U1(4*i + 2) + y1*h_2;

        sin_z = sin(kz*(z+h_2))/rho;
        sin_x = sin(kx*x);
        sinh_y = sinh(ky*y);
        cos_x = cos(kx*x);
        cosh_y = cosh(ky*y);

        Ax = 1./kz * cos_x*cosh_y*sin_z;
        Ay = kx/(ky*kz) * sin_x*sinh_y*sin_z;
        A1 = kx/kz * sin_x*cosh_y*sin_z;
        A2 = cos_x*sinh_y*sin_z;

        x1 = px - Ax;
        y1 = py - Ay;
        Px1 = -x1*A1 + y1*kx*kx/(ky*kz)*A2;
        Py1 = y1*A1 + x1*ky/kz*A2;

        U1(4*(i+1) + 1) = U1(4*i + 1) + Px1*h;
        U1(4*(i+1) + 3) = U1(4*i + 3) + Py1*h;
        U1(4*(i+1) + 0) = U1(4*i + 0) + x1*h;
        U1(4*(i+1) + 2) = U1(4*i + 2) + y1*h;

        z += h;
        }
    """
    weave.inline(code, ['u',"q"])
    x = u[::4]
    y = u[2::4]
    #print x
    #print y
    px = u[1::4] - 1/kz * np.cos(kx*x)*np.cosh(ky*y)*np.sin(kz*z)/rho
    py = u[3::4] - kx/(ky*kz) * np.sin(kx*x)*np.sinh(ky*y)*np.sin(kz*z)/rho
    return x, px,y, py



if __name__ == "__main__":
    from pylab import *
    y0 = [0.1,0.0001,0.002,0.0002,0,0]
    lu = 0.007
    B0 = 0.66
    gamma = 5000
    c = 299792458
    m0 = 0.510998928*1e+6
    kz = 2*pi/lu
    kx = 2*pi/0.65
    Kx =  B0*lu*c/(m0*2.*pi)
    rho = gamma/Kx/kz
    print "x' = ", 1/rho/kz
    print "x = ", 1/rho/kz/kz
    z = linspace(0,lu*200, num = 10000)
    x, px,y, py, z1, pz = track_und_RK(y0, z, kz, kx ,B0, gamma)
    print "x/px on the end: ", x[-1], px[-1]
    print "y/py on the end: ", y[-1], py[-1]
    figure(1)
    plot(z, x, z1, x)
    figure(2)
    plot(z, py)
    show()