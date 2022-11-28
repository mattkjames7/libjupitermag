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

	/* this is for polarization */
	int nalpha = 1;
	double alpha = 0.0;

	/* model names */
	const char *IntFunc = "jrm33";
	int nExt = 1;
	char **ExtFuncs = new char*[nExt];
	ExtFuncs[0] = new char[8];
	strcpy(ExtFuncs[0],"con2020");

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
	double **x = new double*[n];
	x[0] = new double[1000];
	double **y = new double*[n];
	y[0] = new double[1000];
	double **z = new double*[n];
	z[0] = new double[1000];

	double **Bx = new double*[n];
	Bx[0] = new double[1000];
	double **By = new double*[n];
	By[0] = new double[1000];
	double **Bz = new double*[n];
	Bz[0] = new double[1000];

	double **R = new double*[n];
	R[0] = new double[1000];
	double **S = new double*[n];
	S[0] = new double[1000];
	double **Rnorm = new double*[n];
	Rnorm[0] = new double[1000];
	double **FP = new double*[n];
	FP[0] = new double[7];

	int nalpha = 2;
	double alpha[] = {0.0,90.0};
	double halpha[nalpha*n*MaxLen];

	/* call the trace wrapper */
	TraceField(n,x0,y0,z0,IntFunc,nExt,ExtFunc,MaxLen,MaxStep,InitStep,minStep,ErrMax,
			Delta,Verbose,TraceDir,nstep,x,y,z,Bx,By,Bz,R,S,Rnorm,FP,nalpha,alpha,halpha);

	printf("C Test Complete\n")

}