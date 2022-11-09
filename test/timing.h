#ifndef __TIMING_H__
#define __TIMING_H__
#include <stdio.h>
#include <stdlib.h>
#include "../include/jupitermag.h"
#include <ctime>

void ReadDataFile(int *n, double *r, double *t, double *p);
void FreeData(double *r, double *t, double *p);


double TimeCon2020Vectorrtp(int n, double *r, double *t, double *p, const char *eqtype);
double TimeCon2020Vectorxyz(int n, double *x, double *y, double *z, const char *eqtype);
double TimeCon2020Scalarrtp(int n, double *r, double *t, double *p, const char *eqtype);
double TimeCon2020Scalarxyz(int n, double *x, double *y, double *z, const char *eqtype);
double TimeInternalVectorrtp(int n, double *r, double *t, double *p, const char *modelname);
double TimeInternalVectorxyz(int n, double *x, double *y, double *z, const char *modelname);
double TimeInternalScalarrtp(int n, double *r, double *t, double *p, const char *modelname);
double TimeInternalScalarxyz(int n, double *x, double *y, double *z, const char *modelname);

void rtptoxyz(int n, double *r, double *t, double *p,
				double *x, double *y, double *z);
				
#endif
