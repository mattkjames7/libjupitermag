#ifndef __FOOTPRINT_H__
#define __FOOTPRINT_H__
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#include "coordconv.h"

const double deg2rad = M_PI/180.0;
const double rad2deg = 180.0/M_PI;

typedef struct FPstr {
	/* SIII North */
	double xn3;
	double yn3;
	double zn3;
	/* SIII South */
	double xs3;
	double ys3;
	double zs3;
	/* MAG North */
	double xnm;
	double ynm;
	double znm;
	/* MAG South */
	double xsm;
	double ysm;
	double zsm;

	/* SIII North */
	double lonn;
	double latn;
	/* MAG North */
	double mlonn;
	double mlatn;
	/* SIII South */
	double lons;
	double lats;
	/* MAG South */
	double mlons;
	double mlats;

} FPstr;

typedef struct EqFPstr {
	/* SIII coords */
	double x3;
	double y3;
	double z3;

	/* MAG coords */
	double xm;
	double ym;
	double zm;

	/*Equatorial footprint*/
	double lshell;
	double mlone;
	double fllen;
} EqFPstr;


void FillFPOutputArray(int n, FPstr *fpi, FPstr *fps, EqFPstr *fpe, double **FP);

void _fillfp(FPstr *fp, double xt, double xp);


void _getbegfp(	int n, double *x, double *y, double *z,
				bool begfp, bool endfp, int indmxr,
				double a, double b,
				double *xfp, double *yfp, double *zfp);

void _getendfp(	int n, double *x, double *y, double *z,
				bool begfp, bool endfp, int indmxr,
				double a, double b,
				double *xfp, double *yfp, double *zfp);

bool _posisfp(double x, double y, double z, double a, double b);

int _maxR(double n, double *x, double *y, double *z);

void _nsends(int n, double *x, double *y, double *z,
			double xt, double xp,
			bool begfp, bool endfp,
			int *begns, int *endns, int *indmxr);

double _fllen(int n, double *x, double *y, double *z) ;

void calculateEquatorialFootprints(int n, double *x, double *y, double *z,
						double xt, double xp, EqFPstr *efp);

void calculateFootprints(int n, double *x, double *y, double *z,
							double a, double b, double xt, double xp,
							FPstr *fp);

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
