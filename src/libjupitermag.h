#ifndef __LIBJUPITERMAG_H__
#define __LIBJUPITERMAG_H__
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "trace.h"
#include <string.h>
#include "../lib/libcon2020/src/libcon2020.h"
#include "../lib/libinternalfield/include/internalfield.h"



extern "C" {
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
