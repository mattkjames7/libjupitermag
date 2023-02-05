#ifndef __COORDCONV_H__
#define __COORDCONV_H__
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

void MagtoSIII(	double xm, double ym, double zm, 
				double xt, double xp,
				double *x3, double *y3, double *z3);

void SIIItoMag(	double x3, double y3, double z3, 
				double xt, double xp,
				double *xm, double *ym, double *zm);

#endif
