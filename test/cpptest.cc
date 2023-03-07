#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <jupitermag.h>

int main() {
	
	printf("C++ Test...\n");
	int n = 1;
	double x0 = 5.0;
	double y0 = 0.0;
	double z0 = 0.0;
	int nalpha = 1;
	double alpha = 0.0;
	
	printf("Create field function vectors\n");
	/* test function for debugging the trace */
	std::vector<FieldFuncPtr> Funcs;

	/* internal model */
	Funcs.push_back(jrm09Field);

	/* external model */
	Funcs.push_back(Con2020Field);

	/* initialise the trace object */
	printf("Create Trace object\n");
	Trace T(Funcs);

	/* add the starting posiutions fo the traces */
	printf("Add starting position\n");
	T.InputPos(n,&x0,&y0,&z0);

	/* configure the trace parameters */
	printf("Set the trace parameters \n");
	T.SetTraceCFG();

	/* set up the alpha calculation */
	printf("Initialize alpha\n");
	T.SetAlpha(nalpha,&alpha);
	

	/* Trace */
	printf("Trace\n");
	T.TraceField();

	/* trace distance, footprints, Rnorm */
	printf("Footprints etc...\n");
	T.CalculateTraceDist();
	T.CalculateTraceFP();
	T.CalculateTraceRnorm();
	

	/* halpha */
	printf("H_alpha\n");
	T.CalculateHalpha();
	
	printf("Trace s and h_alpha\n");

	double *FP = new double[49];
	T.GetTraceFootprints(&FP);

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

	delete[] FP;
	printf("C++ Test Complete\n");

	

	
}
