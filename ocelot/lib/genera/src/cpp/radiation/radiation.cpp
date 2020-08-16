
#include "definition.h"
//#include "integral.h"
//#include "radiation.h"
//#include "spline.h"
#include <math.h>
#include <omp.h>


int radiation(const Particle particle, const Motion *motion, Screen *screen)
{
    if(motion->bRough>0)
    {
        return 0;
    }
    int ret1 = 1;
    double hc = 1.239841874330e-3; // mm
    double pi = 3.141592653589793;
    double k2q3 = 1.1547005383792517;//  = 2./sqrt(3)
    double Q = 0.5866740802042227; // (mm*T)^-1
    double lperiod = screen->lperiod;
    int nperiods = screen->nperiods;
    int status = screen->status;
    double Zshift = 0.;
    int ncenter = nperiods/2;
    if (status == 13)
    Zshift = lperiod*(ncenter);
    
    
    
    double gamma = particle.gamma;
    double gamma2 = particle.gamma*particle.gamma;
	
    int Nmotion = (motion->size + 1)/3;
    double h2 = (motion->Z[motion->size-1]-motion->Z[0])/2./(Nmotion-1);
    
    double w[] = {0.5555555555555556*h2, 0.8888888888888889*h2, 0.5555555555555556*h2};
    
    int xpoint = screen->xNstep;
    int ypoint = screen->yNstep;
	int crtl = xpoint*ypoint*screen->eNstep;
    
    double Zin = motion->Z[0]+Zshift;
    double Xin = motion->X[0];
    double Yin = motion->Y[0];
    
    int i = 0;
    int n = 0;
    int p = 0;
    
    double XX, YY, ZZ;
    double BetX, BetY;
    double IbetX2, IbetY2;
    double Bx, By;
    double LenPntrZ;
    double radConstAdd;
	
	double LenPntrConst = screen->Distance - Zin; // I have to pay attention to this
    double *reX = new double [crtl];
    double *imX = new double [crtl];
    double *reY = new double [crtl];
    double *imY = new double [crtl];
    //int  nthreats =

    if(xpoint*ypoint*screen->eNstep<600)
    {
        omp_set_num_threads(1);
        omp_set_dynamic(0);
    }
    else
    {
        omp_set_num_threads(16);
        omp_set_dynamic(1);
    }

    for(n = 0; n<Nmotion-1; n++)
    {
        for(p = 0; p<3; p++) // Gauss integration
        {
            
            i = n*3 + p + 1;
            XX = motion->X[i];
            YY = motion->Y[i];
            ZZ = motion->Z[i]+Zshift;
            BetX = motion->Xbeta[i];
            BetY = motion->Ybeta[i];
            IbetX2 = motion->XbetaI2[i];
            IbetY2 = motion->YbetaI2[i];
            Bx = motion->Bx[i];
            By = motion->By[i];
            LenPntrZ = screen->Distance - ZZ; // denominator for pointer nx(z) and ny(z)
            radConstAdd = w[p]*Q*k2q3*(screen->Distance -screen->Zstart);//LenPntrConst;
            
            //omp_set_dynamic(1);      // Á‡ÔÂÚËÚ¸ ·Ë·ÎËÓÚÂÍÂ openmp ÏÂÌˇÚ¸ ˜ËÒÎÓ ÔÓÚÓÍÓ‚ ‚Ó ‚ÂÏˇ ËÒÔÓÎÌÂÌËˇ
            //omp_set_num_threads(8); // ÛÒÚ‡ÌÓ‚ËÚ¸ ˜ËÒÎÓ ÔÓÚÓÍÓ‚ ‚ 10
            #pragma omp parallel for
            for(int jx = 0; jx < xpoint; jx++)
            {
                for(int jy = 0; jy < ypoint; jy++)
                {
                    //#pragma omp parallel for
                    for(int je = 0; je < screen->eNstep; je++)
                    {
                        
                        double Xscr = screen->xStart + screen->xStep*jx;
                        double Yscr = screen->yStart + screen->yStep*jy;
                        double Erad = screen->eStart + screen->eStep*je; // energy !!!!!
                        
                        
                        double prX = Xscr - XX; //for pointer nx(z)
                        double prY = Yscr - YY; //for pointer ny(z)
                        double nx = prX/LenPntrZ;
                        double ny = prY/LenPntrZ;
                        double tx = gamma*(nx - BetX);
                        double ty = gamma*(ny - BetY);
                        double tx2 = tx*tx;
                        double ty2 = ty*ty;
                        double tyx = 2.*tx*ty;
                        double denominator = (1. + tx2 + ty2)*(1. + tx2 + ty2);
                        double radConst = radConstAdd/(LenPntrZ*denominator);
                        /*sigma*/         double radX = radConst*(By*(1. - tx2 + ty2) + Bx*tyx - 2.*tx/Q/LenPntrZ);
                        /*pi*/            double radY =-radConst*(Bx*(1. + tx2 - ty2) + By*tyx + 2.*ty/Q/LenPntrZ);
                        
                        double faseConst = pi*Erad/(gamma2*hc);
                        double prXconst = Xscr - Xin;
                        double prYconst = Yscr - Yin;
                        
                        double phaseConstIn = prXconst*prXconst/LenPntrConst + prYconst*prYconst/LenPntrConst;
                        double phaseConstCur = prX*prX/LenPntrZ + prY*prY/LenPntrZ;
                        // string below is for case direct accumulation
                        //double fase = screen->Phase[ypoint*xpoint*je + xpoint*jy + jx] + faseConst*(ZZ - motion->Z[0]  + gamma2*(IbetX2 + IbetY2 + phaseConstCur - phaseConstIn));
                        double fase = faseConst*(ZZ - Zin + gamma2*(IbetX2 + IbetY2 + phaseConstCur - phaseConstIn));
                        double cosf = cos(fase);
                        double sinf = sin(fase);
                        
                        double EreX = radX*cosf; //(cosf *cos(fase0) - sinf*sin(fase0));
                        double EimX = radX*sinf; //(sinf *cos(fase0) + cosf*sin(fase0));
                        double EreY = radY*cosf;
                        double EimY = radY*sinf;
                        
                        screen->ReEx[ypoint*xpoint*je + xpoint*jy + jx] += EreX;
                        screen->ImEx[ypoint*xpoint*je + xpoint*jy + jx] += EimX;
                        screen->ReEy[ypoint*xpoint*je + xpoint*jy + jx] += EreY;
                        screen->ImEy[ypoint*xpoint*je + xpoint*jy + jx] += EimY;
                        
                        if (n == Nmotion-2 && p == 2) //(n == 5000 && p == 2)
                        {
                            int j = n*3 + p + 2;
                            prX = Xscr - motion->X[j]; //for pointer nx(z)
                            prY = Yscr - motion->Y[j]; //for pointer ny(z)
                            IbetX2 = motion->XbetaI2[j];
                            IbetY2 = motion->YbetaI2[j];
                            fase = faseConst*(motion->Z[j] - motion->Z[0]  + gamma2*(IbetX2 + IbetY2 + prX*prX/LenPntrZ + prY*prY/LenPntrZ - phaseConstIn));
                            screen->Phase[ypoint*xpoint*je + xpoint*jy + jx] = fase;
                            if(status == 13)
                            {
                                reX[ypoint*xpoint*je + xpoint*jy + jx] = screen->ReEx[ypoint*xpoint*je + xpoint*jy + jx];
                                imX[ypoint*xpoint*je + xpoint*jy + jx] = screen->ImEx[ypoint*xpoint*je + xpoint*jy + jx];
                                reY[ypoint*xpoint*je + xpoint*jy + jx] = screen->ReEy[ypoint*xpoint*je + xpoint*jy + jx];
                                imY[ypoint*xpoint*je + xpoint*jy + jx] = screen->ImEy[ypoint*xpoint*je + xpoint*jy + jx];
                                //phase_array[ypoint*xpoint*je + xpoint*jy + jx] = screen->Phase[ypoint*xpoint*je + xpoint*jy + jx];
                                screen->ReEx[ypoint*xpoint*je + xpoint*jy + jx] = 0.;
                                screen->ImEx[ypoint*xpoint*je + xpoint*jy + jx] = 0.;
                                screen->ReEy[ypoint*xpoint*je + xpoint*jy + jx] = 0.;
                                screen->ImEy[ypoint*xpoint*je + xpoint*jy + jx] = 0.;
                                screen->Phase[ypoint*xpoint*je + xpoint*jy + jx] = 0.;
                                //phase_end = fase;
                                //ZZ_end_calc = motion->Z[j];
                            }
                            
                        }
                        
                    }
                }
            }
        } // Gauss integration
    }
    
    if(status == 13)
    {
        double IbetX2_end_period = IbetX2;
        double IbetY2_end_period = IbetY2;
        Zin = 0.;
        for(p = 1; p<=nperiods; p++)
        {
            ZZ = lperiod*p;
            IbetX2 = IbetX2_end_period*(p);
            IbetY2 = IbetY2_end_period*(p);
            
            LenPntrZ = screen->Distance - ZZ; // denominator for pointer nx(z) and ny(z)
            //radConstAdd = w[p]*Q*k2q3*(screen->Distance - screen->Zstart);//LenPntrConst;
            
            //omp_set_dynamic(0);      // Á‡ÔÂÚËÚ¸ ·Ë·ÎËÓÚÂÍÂ openmp ÏÂÌˇÚ¸ ˜ËÒÎÓ ÔÓÚÓÍÓ‚ ‚Ó ‚ÂÏˇ ËÒÔÓÎÌÂÌËˇ
            //omp_set_num_threads(8); // ÛÒÚ‡ÌÓ‚ËÚ¸ ˜ËÒÎÓ ÔÓÚÓÍÓ‚ ‚ 10
            #pragma omp parallel for
            for(int jx = 0; jx < xpoint; jx++)
            {
                for(int jy = 0; jy < ypoint; jy++)
                {
                    for(int je = 0; je < screen->eNstep; je++)
                    {
                        
                        double Xscr = screen->xStart + screen->xStep*jx;
                        double Yscr = screen->yStart + screen->yStep*jy;
                        double Erad = screen->eStart + screen->eStep*je; // energy !!!!!
                        
                        
                        //double prX = Xscr; //for pointer nx(z)
                        //double prY = Yscr; //for pointer ny(z)
                        
                        double faseConst = pi*Erad/(gamma2*hc);
                        double prXconst = Xscr - Xin;
                        double prYconst = Yscr - Yin;
                        
                        double phaseConstIn = prXconst*prXconst/LenPntrConst + prYconst*prYconst/LenPntrConst;
                        double phaseConstCur = Xscr*Xscr/LenPntrZ + Yscr*Yscr/LenPntrZ;
                        //string below is for case direct accumulation
                        //double fase = screen->Phase[ypoint*xpoint*je + xpoint*jy + jx] + faseConst*(ZZ - motion->Z[0]  + gamma2*(IbetX2 + IbetY2 + phaseConstCur - phaseConstIn));
                        double fase = faseConst*(ZZ - Zin + gamma2*(IbetX2 + IbetY2 + phaseConstCur - phaseConstIn));
                        double delta_phase = fase - screen->Phase[ypoint*xpoint*je + xpoint*jy + jx];
                        double cosf = cos(fase-delta_phase);
                        double sinf = sin(fase-delta_phase);
                        
                        double normalize = (screen->Distance - Zshift-lperiod)/(LenPntrZ);
                        double EreX = reX[ypoint*xpoint*je + xpoint*jy + jx]*normalize; //(cosf *cos(fase0) - sinf*sin(fase0));
                        double EimX = imX[ypoint*xpoint*je + xpoint*jy + jx]*normalize; //(sinf *cos(fase0) + cosf*sin(fase0));
                        double EreY = reY[ypoint*xpoint*je + xpoint*jy + jx]*normalize;
                        double EimY = imY[ypoint*xpoint*je + xpoint*jy + jx]*normalize;
                        
                        screen->ReEx[ypoint*xpoint*je + xpoint*jy + jx] += EreX*cosf - EimX*sinf;
                        screen->ImEx[ypoint*xpoint*je + xpoint*jy + jx] += EimX*cosf + EreX*sinf;
                        screen->ReEy[ypoint*xpoint*je + xpoint*jy + jx] += EreY*cosf - EimY*sinf;
                        screen->ImEy[ypoint*xpoint*je + xpoint*jy + jx] += EimY*cosf + EreY*sinf;
                        screen->Phase[ypoint*xpoint*je + xpoint*jy + jx] = fase;
                        //phase_array[ypoint*xpoint*je + xpoint*jy + jx] = fase;
                    }//je
                }//jy
            }//jx
        } // nperiod
    } // if status
    delete []reX;
    delete []imX;
    delete []reY;
    delete []imY;
    return ret1;
}




void sum_screens(Screen *screen, Screen screen_up)
{
    //Screen screen; 
    int size = screen->eNstep*screen->xNstep*screen->yNstep;

    for(int i = 0; i<size; i++)
    {
        double sinfa = sin(screen->Phase[i]);
        double cosfa = cos(screen->Phase[i]);
        screen->ReEx[i] += screen_up.ReEx[i]*cosfa - screen_up.ImEx[i]*sinfa;
        screen->ImEx[i] += screen_up.ImEx[i]*cosfa + screen_up.ReEx[i]*sinfa;
        screen->ReEy[i] += screen_up.ReEy[i]*cosfa - screen_up.ImEy[i]*sinfa;
        screen->ImEy[i] += screen_up.ImEy[i]*cosfa + screen_up.ReEy[i]*sinfa;
        screen->Phase[i] += screen_up.Phase[i];
    }
}




//__declspec(dllexport) void trajectory(const double *arZ, const double *arFz, const int size)
extern "C"  int emsolver(const int numberElem, const int *stepMotion, double **arMotion, const double gamma, const double *scrPrm, double *arScreen)
{
    Particle particle;


    particle.X = 0.;
    particle.Y = 0.;
    particle.Z = 0.;
    particle.betaX = 0.;
    particle.betaY = 0.;
    particle.gamma = gamma;
    particle.mass = 1.;
    particle.charge = -1;
    particle.curr = 1.;

    Screen screen;
    screen.eNstep = scrPrm[0];
    screen.eStart = scrPrm[1];
    screen.eStep = scrPrm[2];

    screen.xNstep = scrPrm[3];
    screen.xStart = scrPrm[4];
    screen.xStep = scrPrm[5];

    screen.yNstep = scrPrm[6];
    screen.yStart = scrPrm[7];
    screen.yStep = scrPrm[8];

    screen.Distance = scrPrm[9];
    screen.Zstart = scrPrm[10];
    screen.lperiod = scrPrm[11];
    screen.nperiods = scrPrm[12];
    screen.status = scrPrm[13];
    //screen.First = scrPrm[11];



    int sizeScr = screen.eNstep*screen.xNstep*screen.yNstep;

    screen.ReEx = &arScreen[0];
    screen.ImEx = &arScreen[sizeScr];
    screen.ReEy = &arScreen[2*sizeScr];
    screen.ImEy = &arScreen[3*sizeScr];
    screen.Phase = &arScreen[4*sizeScr];

    // Screen UP
    Screen screen_up;

    screen_up.ReEx = new double [sizeScr];
    screen_up.ImEx = new double [sizeScr];
    screen_up.ReEy = new double [sizeScr];
    screen_up.ImEy = new double [sizeScr];
    screen_up.Phase = new double [sizeScr];
    screen_up.eNstep = screen.eNstep;
    screen_up.eStart = screen.eStart;
    screen_up.eStep = screen.eStep;

    screen_up.xNstep = screen.xNstep;
    screen_up.xStart = screen.xStart;
    screen_up.xStep = screen.xStep;

    screen_up.yNstep = screen.yNstep;
    screen_up.yStart = screen.yStart;
    screen_up.yStep = screen.yStep;

    screen_up.Distance = screen.Distance;
    screen_up.Zstart = screen.Zstart;
    screen_up.lperiod = scrPrm[11];
    screen_up.nperiods = scrPrm[12];
    screen_up.status = scrPrm[13];
    //screen_up.First = screen.First;
    // Screen UP


    int ret = 1;

    for (int i = 0; i<numberElem; i++)
    {
        for(int n = 0; n<sizeScr; n++)
        {
            screen_up.ReEx[n] = 0.;
            screen_up.ImEx[n] = 0.;
            screen_up.ReEy[n] = 0.;
            screen_up.ImEy[n] = 0.;
            screen_up.Phase[n] = 0.;
        }

        Motion motion;
        motion.bRough = 0;
        int size = stepMotion[i];
        motion.size = size;
        motion.X = &arMotion[i][0];
        motion.Y = &arMotion[i][size];
        motion.Z = &arMotion[i][2*size];
        motion.Xbeta = &arMotion[i][3*size];
        motion.Ybeta = &arMotion[i][4*size];
        motion.Zbeta = &arMotion[i][5*size];
        motion.Bx = &arMotion[i][6*size];
        motion.By = &arMotion[i][7*size];
        motion.Bz = &arMotion[i][8*size];
        motion.XbetaI2 = &arMotion[i][9*size];
        motion.YbetaI2 = &arMotion[i][10*size];

        ret = radiation(particle, &motion, &screen_up);
        sum_screens(&screen, screen_up);

    }
    delete [] screen_up.ReEx;
    delete [] screen_up.ImEx;
    delete [] screen_up.ReEy;
    delete [] screen_up.ImEy;
    delete [] screen_up.Phase;

    return ret;
}
