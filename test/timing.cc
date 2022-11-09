#include "timing.h"

void ReadDataFile(int *n, double *r, double *t, double *p) {

	/* open the file and get the number of elements */
	FILE *f = open('data.bin','r')
	fread(n,1,sizeof(int),f)
	printf("Found %d elements\n",*n);

	/* allocate temporary arrays*/
	float *rf = new float[*n];
	float *tf = new float[*n];
	float *pf = new float[*n];

	/* read them in */
	printf("Reading r\n");
	fread(rf,*n,sizeof(float),f);
	printf("Reading t\n");
	fread(tf,*n,sizeof(float),f);
	printf("Reading p\n");
	fread(pf,*n,sizeof(float),f);
	
	/* close the file */
	fclose(f);

	/* convert to doubles */
	r = new double[*n];
	t = new double[*n];
	p = new double[*n];
	int i;
	for (i=0;i<*n;i++) {
		r[i] = (double) rf[i];
		t[i] = (double) tf[i];
		p[i] = (double) pf[i];
	}

	/* delete the 32 bit floats */
	delete [] rf;
	delete [] tf;
	delete [] pf;

	printf("Done\n");

}

void FreeData(double *r, double *t, double *p) {
	delete [] r;
	delete [] t;
	delete [] p;
}

double TimeCon2020Vectorrtp(int n, double *r, double *t, double *p, const char *eqtype) {

	/* create empty arrays to store B vectors in */
	double *Br = new double[n];
	double *Bt = new double[n];
	double *Bp = new double[n];

	/* create a model object */
	Con2020 model;

	/* initialize it */
	model.SetEqType(eqtype);
	model.SetSmooth(false);
	model.SetCartIn(false);
	model.SetCartOut(false);

	/* time it */
	double t0, t1;
	t0 = clock();
	model.Field(n,r,t,p,Br,Bt,Bp);
	t1 = clock();

	/* delete the allocated field vector arrays */
	FreeData(Br,Bt,Bp);

	return (t1 - t0)/CLOCKS_PER_SEC;

}

double TimeCon2020Vectorxyz(int n, double *x, double *y, double *z, const char *eqtype) {

	/* create empty arrays to store B vectors in */
	double *Bx = new double[n];
	double *By = new double[n];
	double *Bz = new double[n];

	/* create a model object */
	Con2020 model;

	/* initialize it */
	model.SetEqType(eqtype);
	model.SetSmooth(false);
	model.SetCartIn(true);
	model.SetCartOut(true);

	/* time it */
	double t0, t1;
	t0 = clock();
	model.Field(n,x,y,z,Bx,By,Bz);
	t1 = clock();

	/* delete the allocated field vector arrays */
	FreeData(Bx,By,Bz);

	return (t1 - t0)/CLOCKS_PER_SEC;

}

double TimeCon2020Scalarrtp(int n, double *r, double *t, double *p, const char *eqtype) {

	/* create empty variables for the field vectors */
	double Br;
	double Bt;
	double Bp;

	/* create a model object */
	Con2020 model;

	/* initialize it */
	model.SetEqType(eqtype);
	model.SetSmooth(false);
	model.SetCartIn(false);
	model.SetCartOut(false);

	/* time it */
	double t0, t1;
	t0 = clock();
	int i;
	for (i=0;i<n;i++) {
		model.Field(r[i],t[i],p[i],&Br,&Bt,&Bp);
	}
	t1 = clock();

	return (t1 - t0)/CLOCKS_PER_SEC;

}

double TimeCon2020Scalarxyz(int n, double *x, double *y, double *z, const char *eqtype) {

	/* create empty variables for the field vectors */
	double Bx;
	double By;
	double Bz;

	/* create a model object */
	Con2020 model;

	/* initialize it */
	model.SetEqType(eqtype);
	model.SetSmooth(false);
	model.SetCartIn(true);
	model.SetCartOut(true);

	/* time it */
	double t0, t1;
	t0 = clock();
	int i;
	for (i=0;i<n;i++) {
		model.Field(x[i],y[i],z[i],&Bx,&By,&Bz);
	}
	t1 = clock();

	return (t1 - t0)/CLOCKS_PER_SEC;

}



double TimeInternalVectorrtp(int n, double *r, double *t, double *p, const char *modelname) {

	/* create empty arrays to store B vectors in */
	double *Br = new double[n];
	double *Bt = new double[n];
	double *Bp = new double[n];

	/* create a model object */
	Internal model = Internal(modelname);

	/* time it */
	double t0, t1;
	t0 = clock();
	model.Field(n,r,t,p,Br,Bt,Bp);
	t1 = clock();

	/* delete the allocated field vector arrays */
	FreeData(Br,Bt,Bp);

	return (t1 - t0)/CLOCKS_PER_SEC;

}




double TimeInternalVectorxyz(int n, double *x, double *y, double *z, const char *modelname) {

	/* create empty arrays to store B vectors in */
	double *Bx = new double[n];
	double *By = new double[n];
	double *Bz = new double[n];

	/* create a model object */
	Internal model = Internal(modelname);

	/* time it */
	double t0, t1;
	t0 = clock();
	model.Field(n,x,y,z,Bx,By,Bz);
	t1 = clock();

	/* delete the allocated field vector arrays */
	FreeData(Bx,By,Bz);

	return (t1 - t0)/CLOCKS_PER_SEC;

}


double TimeInternalScalarrtp(int n, double *r, double *t, double *p, const char *modelname) {

	/* create empty variables for the field vectors */
	double Br;
	double Bt;
	double Bp;

	/* create a model object */
	Internal model = Internal(modelname);

	/* time it */
	double t0, t1;
	t0 = clock();
	int i;
	for (i=0;i<n;i++) {
		model.Field(r[i],t[i],p[i],&Br,&Bt,&Bp);
	}
	t1 = clock();


	return (t1 - t0)/CLOCKS_PER_SEC;

}


double TimeInternalScalarxyz(int n, double *x, double *y, double *z, const char *modelname) {

	/* create empty variables for the field vectors */
	double Bx;
	double By;
	double Bz;

	/* create a model object */
	Internal model = Internal(modelname);

	/* time it */
	double t0, t1;
	t0 = clock();
	int i;
	for (i=0;i<n;i++) {
		model.Field(x[i],y[i],z[i],&Bx,&By,&Bz);
	}
	t1 = clock();


	return (t1 - t0)/CLOCKS_PER_SEC;

}


void rtptoxyz(int n, double *r, double *t, double *p,
				double *x, double *y, double *z) {
	
	int i;
	for (i=0;i<n;i++) {
		x[i] = r[i]*sin(t[i])*cos(p[i]);
		y[i] = r[i]*sin(t[i])*sin(p[i]);
		z[i] = r[i]*cos(t[i]);
	}
}



int main() {

	/* read in positions to test models with */
	double *r, *t, *p;
	int n;
	ReadDataFile(&n,r,t,p);

	/* get the cartesian positions too */
	rtptoxyz(n,r,t,p,x,y,z);


	double dt;

	/* test the internal model (vector,rtp)*/
	dt = TimeInternalScalarrtp(n,r,t,p,"jrm33");
	printf("Vector, RTP %s: %f s\n","jrm33",dt);

	/* free memory */
	FreeData(r,t,p);
	FreeData(Br,Bt,Bp);
}