#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H
#include "motion.h"
#include "spline.h"
#include "u_magfield.h"
#include "u_field.h"

int rkgauss( const int NstepMotion,
             const Particle particle,
			 Field *field,
             Motion *motion);

#endif // RUNGEKUTTA_H
