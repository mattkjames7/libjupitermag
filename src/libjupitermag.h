#ifndef __LIBJUPITERMAG_H__
#define __LIBJUPITERMAG_H__
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "trace.h"
#include <string.h>
#include "libcon2020.h"
#include "internalfield.h"



extern "C" {
	void JupitermagGetCon2020Params(double *mui, double *irho, double *r0, double *r1,
				double *d, double *xt, double *xp, char *eqtype,
				bool *Edwards, bool *ErrChk, bool *CartIn, bool *CartOut,
				bool *smooth, double *DeltaRho, double *DeltaZ,
				double *g, char *azfunc, double *wO_open, double *wO_om,
				double *thetamm, double *dthetamm, double *thetaoc, double *dthetaoc);

	bool TraceField(int n, double *x0, double *y0, double *z0,
				const char *IntFunc, int nExt, char **ExtFunc,
				int MaxLen, double MaxStep, double InitStep,
				double MinStep, double ErrMax, double Delta,
				bool Verbose, int TraceDir,
				double as, double bs, double ai, double bi,
				int *nstep,
				double **x, double **y, double **z,
				double **Bx, double **By, double **Bz,
				double **R, double **S, double **Rnorm, 
				int **traceRegion, double **FP,
				int nalpha, double *alpha, double *halpha);
}
#endif
