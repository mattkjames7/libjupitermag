#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <jupitermagc.h>

int main () {
	
	printf("C Test...\n");

	/* start of the trace */
	int n = 1;
	double x0 = 5.0;
	double y0 = 0.0;
	double z0 = 0.0;


	/* model names */
	const char *IntFunc = "jrm33";
	int nExt = 1;
	char *ExtFuncs = (char*) malloc(nExt*sizeof(char));
	strcpy(ExtFuncs,"con2020");

	/* trace configuration */
	int MaxLen = 1000;
	double MaxStep = 1.0;
	double InitStep = 0.5;
	double MinStep = 0.001;
	double ErrMax = 0.0001;
	double Delta = 0.05;
	bool Verbose = true;
	int TraceDir = 0;

	/* trace output */
	int nstep[n];
	double *x = (double *) malloc(MaxLen*sizeof(double));
	double *y = (double *) malloc(MaxLen*sizeof(double));
	double *z = (double *) malloc(MaxLen*sizeof(double));

	double *Bx = (double *) malloc(MaxLen*sizeof(double));
	double *By = (double *) malloc(MaxLen*sizeof(double));
	double *Bz = (double *) malloc(MaxLen*sizeof(double));

	double *R = (double *) malloc(MaxLen*sizeof(double));
	double *S = (double *) malloc(MaxLen*sizeof(double));
	double *Rnorm = (double *) malloc(MaxLen*sizeof(double));
	double *FP = (double *) malloc(MaxLen*sizeof(double));

	/* this is for polarization */
	int nalpha = 2;
	double alpha[] = {0.0,90.0};
	double halpha[nalpha*n*MaxLen];

	/* call the trace wrapper */
	TraceField(n,&x0,&y0,&z0,IntFunc,nExt,&ExtFuncs,MaxLen,MaxStep,InitStep,MinStep,ErrMax,
			Delta,Verbose,TraceDir,nstep,&x,&y,&z,&Bx,&By,&Bz,&R,&S,&Rnorm,&FP,nalpha,alpha,halpha);

	printf("C Test Complete\n");

}