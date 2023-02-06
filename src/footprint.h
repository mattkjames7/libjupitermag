#ifndef __FOOTPRINT_H__
#define __FOOTPRINT_H__
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "footprint.h"
#include "coordconv.h"


void footprints(int n, double *x, double *y, double *z, 
				double a, double b, double xt, double xp,
				double *xfn, double *yfn, double *zfn,
				double *xfs, double *yfs, double *zfs);

double planetRadius(double x, double y, double z, double a, double b);

double thetaPlanetRadius(double theta, double a, double b);

bool isCrossing(double x0, double y0, double z0,
				double x1, double y1, double z1,
				double a, double b);

void interpCrossing(double x0, double y0, double z0,
					double x1, double y1, double z1,
					double a, double b,
					double *xfp, double *yfp, double *zfp);				

void findFootprint(	double *x, double *y, double *z,
					int starti, int endi, 
					double a, double b,
					double *xfp, double *yfp, double *zfp);

void eqfootprints(	int n, double *x, double *y, double *z,
					double *xfe, double *yfe, double *zfe, 
					double *L, double *Lon);

#endif
