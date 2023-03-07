#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <jupitermag.h>

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
	char *ExtFuncs = (char*) malloc(32*sizeof(char));
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
	double *FP = (double *) malloc(49*sizeof(double));

	/* this is for polarization */
	int nalpha = 2;
	double alpha[] = {0.0,90.0};
	double halpha[nalpha*n*MaxLen];

	/* call the trace wrapper */
	TraceField(n,&x0,&y0,&z0,IntFunc,nExt,&ExtFuncs,MaxLen,MaxStep,InitStep,MinStep,ErrMax,
			Delta,Verbose,TraceDir,nstep,&x,&y,&z,&Bx,&By,&Bz,&R,&S,&Rnorm,&FP,nalpha,alpha,halpha);
	printf("Trace Done...\n");
	printf("**** North Ionospheric FP ****\n");
	printf("lat: %f lon %f\n",FP[31],FP[30]);
	printf("xyz: [ %f %f %f ]\n",FP[0],FP[1],FP[2]);
	printf("mlat: %f mlon %f\n",FP[33],FP[32]);
	printf("xyz: [ %f %f %f ]\n",FP[6],FP[7],FP[8]);


	printf("**** South Ionospheric FP ****\n");
	printf("lat: %f lon %f\n",FP[35],FP[34]);
	printf("xyz: [ %f %f %f ]\n",FP[3],FP[4],FP[5]);
	printf("mlat: %f mlon %f\n",FP[37],FP[36]);
	printf("xyz: [ %f %f %f ]\n",FP[9],FP[10],FP[11]);

	printf("**** North Surface FP ****\n");
	printf("lat: %f lon %f\n",FP[39],FP[38]);
	printf("xyz: [ %f %f %f ]\n",FP[12],FP[13],FP[14]);
	printf("mlat: %f mlon %f\n",FP[41],FP[40]);
	printf("xyz: [ %f %f %f ]\n",FP[18],FP[19],FP[20]);


	printf("**** South Surface FP ****\n");
	printf("lat: %f lon %f\n",FP[43],FP[42]);
	printf("xyz: [ %f %f %f ]\n",FP[15],FP[16],FP[17]);
	printf("mlat: %f mlon %f\n",FP[45],FP[44]);
	printf("xyz: [ %f %f %f ]\n",FP[21],FP[22],FP[23]);

	printf("**** Equatorial FP ****\n");
	printf("mlon: %f lshell %f\n",FP[47],FP[46]);
	printf("xyz: [ %f %f %f ]\n",FP[27],FP[28],FP[29]);


	printf("C Test Complete\n");

	free(ExtFuncs);
	free(x);
	free(y);
	free(z);
	free(Bx);
	free(By);
	free(Bz);
	free(R);
	free(S);
	free(Rnorm);
	free(FP);

}