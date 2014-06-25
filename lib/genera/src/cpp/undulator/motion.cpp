#include "motion.h"
struct Motion;
void initMotion(Motion *motion, const int bRough, const int NstepMotion,int Nsuper, int startPosition, double *aMotion)
{
	motion->bRough = bRough; // change if you want to calculate trajectory for gauss integration 
	motion->size = NstepMotion; 
	int lenArr = Nsuper*NstepMotion - (Nsuper - 1)*(2-bRough); // because I taking into accaunt opportunity use period of undulator

	motion->X = &aMotion[startPosition];
	motion->Y = &aMotion[lenArr + startPosition];
	motion->Z = &aMotion[2*lenArr + startPosition];
	motion->Xbeta = &aMotion[3*lenArr + startPosition];
	motion->Ybeta = &aMotion[4*lenArr + startPosition];
	motion->Zbeta = &aMotion[5*lenArr + startPosition];
	motion->Bx = &aMotion[6*lenArr + startPosition];
	motion->By = &aMotion[7*lenArr + startPosition];
	motion->Bz = &aMotion[8*lenArr + startPosition];

	motion->XbetaI2 = &aMotion[9*lenArr + startPosition];
	motion->YbetaI2 = &aMotion[10*lenArr + startPosition];
}