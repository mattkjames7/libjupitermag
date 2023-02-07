#include "trace.h"

Trace::Trace(std::vector<FieldFuncPtr> Funcs) {
	
	/* set the field fucntion pointers vector */
	Funcs_ = Funcs;
	nf_ = Funcs.size();

	/* initialize all of the boolean parameters */
	inputPos_ = false;
	tracedField_ = false;
	allocTrace_ = false;
	hasFootprints_ = false;
	allocFootprints_ = false;
	hasDist_ = false;
	allocDist_ = false;
	hasRnorm_= false;
	allocRnorm_ = false;
	hasHalpha_= false;
	allocHalpha_ = false;
	allocHalpha3D_ = false;
	allocEqFP_ = false;
	allocAlpha_ = false;
	xt_ = 10.25*deg2rad;
	xp_ = 163.62*deg2rad;
	
	
	/* default trace parameters */
	SetTraceCFG();
	
	SetTraceBoundDefaults();

	/* this is a hack - used to end trace at lower altitudes, while
	 * retaining Jupiter's oblateness */
	RMultiplier_ = 1.0;
		
}

Trace::~Trace() {

	/* check for each allocated variable and delete it*/
	int i, j, k = 0;
	
	/* starting positions */
	if (inputPos_) {
		delete[] x0_;
		delete[] y0_;
		delete[] z0_;
	}

	/* traces */
	if (allocTrace_) {
		for (i=0;i<n_;i++) {
			delete[] x_[i];
			delete[] y_[i];
			delete[] z_[i];
			delete[] bx_[i];
			delete[] by_[i];
			delete[] bz_[i];
			delete[] R_[i];
		}
		delete[] x_;
		delete[] y_;
		delete[] z_;
		delete[] bx_;
		delete[] by_;
		delete[] bz_;
		delete[] R_;
		delete[] nstep_;
	}

	/* field line footprints */
	if (allocFootprints_) {
		for (i=0;i<n_;i++) {
			delete[] FP_[i];
		}
		delete[] FP_;
	}

	/* field line distance */
	if (allocDist_) {
		for (i=0;i<n_;i++) {
			delete[] S_[i];
		}
		delete[] S_;
	}

	/* r norm distance */
	if (allocRnorm_) {
		for (i=0;i<n_;i++) {
			delete[] Rnorm_[i];
		}
		delete[] Rnorm_;
	}

	/* h alpha*/
	if (allocAlpha_) {
		delete[] alpha0_;
		delete[] alpha1_;
	}
	if (allocHalpha_) {
		delete[] Halpha_;
	}
	if (allocHalpha3D_) {
		for (i=0;i<n_;i++) {
			for (j=0;j<nalpha_;j++) {
				delete[] Halpha3D_[i][j];
			}
			delete[] Halpha3D_[i];
		}
		delete[] Halpha3D_;
	}

	/* footprint/endpoints */
	if (allocEqFP_) {
		delete[] xfn3_;
		delete[] yfn3_;
		delete[] zfn3_;
		delete[] xfs3_;
		delete[] yfs3_;
		delete[] zfs3_;
		delete[] xin3_;
		delete[] yin3_;
		delete[] zin3_;
		delete[] xis3_;
		delete[] yis3_;
		delete[] zis3_;
		delete[] xfe3_;
		delete[] yfe3_;
		delete[] zfe3_;

		delete[] xfnm_;
		delete[] yfnm_;
		delete[] zfnm_;
		delete[] xfsm_;
		delete[] yfsm_;
		delete[] zfsm_;
		delete[] xinm_;
		delete[] yinm_;
		delete[] zinm_;
		delete[] xism_;
		delete[] yism_;
		delete[] zism_;
		delete[] xfem_;
		delete[] yfem_;
		delete[] zfem_;
	}
}

void Trace::InputPos(	int n, double *x, double *y, double *z) {

	/* check that we have not already done this */
	if (inputPos_) {
		printf("Input positions already set, ignoring...\n");
		return;
	}

	/* allocate the memory to store the input coords */
	n_ = n;
	x0_ = new double[n];
	y0_ = new double[n];
	z0_ = new double[n];
	int i;
	for (i=0;i<n_;i++) {
		x0_[i] = x[i];
		y0_[i] = y[i];
		z0_[i] = z[i];
	}
	
	/* set the flag so we know to delete  again once the object is deleted */
	inputPos_ = true;
						
}		

void Trace::SetTraceCFG(int MaxLen, double MaxStep, double InitStep,
						double MinStep, double ErrMax, double Delta,
						bool Verbose, int TraceDir) {
	
	/* set all of the params */					
	MaxLen_ = MaxLen;
	MaxStep_ = MaxStep;
	MinStep_ = MinStep;
	InitStep_ = InitStep;
	Verbose_ = Verbose;
	TraceDir_ = TraceDir;
	ErrMax_ = ErrMax;
	Delta_ = Delta;
	MaxR_ = 1000.0;
}

void Trace::SetTraceCFG() {
	
	/* set default params */					
	MaxLen_ = 1000;
	MaxStep_ = 1.0;
	InitStep_ = 0.5;
	MinStep_ = 0.001;
	Verbose_ = false;
	TraceDir_ = 0;
	ErrMax_ = 0.0001;
	Delta_ = 0.05;
	MaxR_ = 1000.0;
}

void Trace::SetAlpha(int nalpha, double *alpha) {
	
	/*NOTE: for each alpha, there will be two traces - one for the 
	 * supplied value and one for alpha + 180 */
	/* set the alpha pointer */
	nalpha_ = nalpha;
	//alpha_ = alpha;
	if (nalpha > 0) {
		alpha0_ = new double[nalpha_];
		alpha1_ = new double[nalpha_];
		allocAlpha_ = true;
		double dtor = M_PI/180.0;
		int i;
		for (i=0;i<nalpha;i++) {
			alpha0_[i] = alpha[i]*dtor;
			alpha1_[i] = fmod(alpha[i]*dtor + M_PI,2*M_PI);
		}
	}
}


Trace Trace::TracePosition(int i, double x, double y, double z) {
	/* return a new trace object at the supplied position using the
	 * parameters at time i */
	Trace T(Funcs_);
	
	/* input position and time - I am pretty certain that the midpoints
	 * of the field lines are stored in SM coords */
	T.InputPos(1,&x,&y,&z);
	
	/* set the model up */
	T.SetTraceCFG(MaxLen_,MaxStep_,InitStep_,MinStep_,ErrMax_,Delta_,false,0);
	T.RMultiplier_ = 0.8;
	
	/* run the GSM trace */
	T.TraceField();
	
	/* calculate S*/
	T.CalculateTraceDist();
	
	return T;
	
}

void Trace::_CalculateTraceHalpha(	int i, int j, double *halpha) {

	/* some variables needed */
	double xe0,ye0,ze0,xe1,ye1,ze1;

	/* get the trace starting points first */
	_CalculateHalphaStartPoints(i,j,&xe0,&ye0,&ze0,&xe1,&ye1,&ze1);

	/* do two traces */
	Trace T0 = TracePosition(i,xe0,ye0,ze0);
	Trace T1 = TracePosition(i,xe1,ye1,ze1);
	
	/* get the closest points to each step of the original trace*/
	double *xc0 = new double[nstep_[i]];
	double *yc0 = new double[nstep_[i]];
	double *zc0 = new double[nstep_[i]];
	double *xc1 = new double[nstep_[i]];
	double *yc1 = new double[nstep_[i]];
	double *zc1 = new double[nstep_[i]];

	interptraceClosestPos(	nstep_[i],x_[i],y_[i],z_[i],
						bx_[i],by_[i],bz_[i],
						T0.nstep_[0],T0.x_[0],T0.y_[0],T0.z_[0],T0.S_[0],
						T1.nstep_[0],T1.x_[0],T1.y_[0],T1.z_[0],T1.S_[0],
						xc0,yc0,zc0,xc1,yc1,zc1);

	/* calculate distances and then halpha */
	double d, dx, dy, dz, h0, h1;
	int k;
	for (k=0;k<nstep_[i];k++) {
		dx = x_[i][k] - xc0[k];
		dy = y_[i][k] - yc0[k];
		dz = z_[i][k] - zc0[k];
		d = sqrt(dx*dx + dy*dy + dz*dz);
		h0 = d/Delta_;
		
		dx = x_[i][k] - xc1[k];
		dy = y_[i][k] - yc1[k];
		dz = z_[i][k] - zc1[k];
		d = sqrt(dx*dx + dy*dy + dz*dz);
		h1 = d/Delta_;
		
		halpha[k] = 0.5*(h0 + h1);
		//printf("%d %5.3f %5.3f %5.3f\n",k,h0,h1,halpha[k]);
	}
	/* free up memory */
	delete[] xc0;
	delete[] yc0;
	delete[] zc0;
	delete[] xc1;
	delete[] yc1;
	delete[] zc1;
}

void Trace::_CalculateHalpha() {

	/* loop through each trace and alpha combination */
	int i, j, k, I, J;
	for (i=0;i<n_;i++) {
		I = i*(nalpha_*MaxLen_);
		for (j=0;j<nalpha_;j++) {
			J = j*MaxLen_;
			_CalculateTraceHalpha(i,j,Halpha3D_[i][j]);
			for (k=0;k<MaxLen_;k++) {
				Halpha_[I + J + k] = Halpha3D_[i][j][k];
			}
		}
	}
}
		
		
bool Trace::_CheckHalpha() {
	
	if (!allocAlpha_) {
		printf("Run the 'SetAlpha()' function prior to calculating h_alpha\n");
		return false;
	}
	
	if (nalpha_ <= 0) {
		printf("1 or more values of alpha must be provided to calculate h_alpha\n");
		return false;
	}
	
	return true;
	
}
		
void Trace::CalculateHalpha() {

	if (!_CheckHalpha()) {
		return;
	}
	
	/* allocate both 1D and 3D arrays */
	Halpha_ = new double[n_*nalpha_*MaxLen_];
	Halpha3D_ = new double**[n_];
	int i, j;
	for (i=0;i<n_;i++) {
		Halpha3D_[i] = new double*[nalpha_];
		for (j=0;j<nalpha_;j++) {
			Halpha3D_[i][j] = new double[MaxLen_];
		}
	}
	allocHalpha_ = true;
	allocHalpha3D_ = true;

	_CalculateHalpha();
}

void Trace::CalculateHalpha(double *halpha) {

	if (!_CheckHalpha()) {
		return;
	}
	
	/* allocate 3D array and use pointer for 1D */
	Halpha_ = halpha;
	Halpha3D_ = new double**[n_];
	int i, j;
	for (i=0;i<n_;i++) {
		Halpha3D_[i] = new double*[nalpha_];
		for (j=0;j<nalpha_;j++) {
			Halpha3D_[i][j] = new double[MaxLen_];
		}
	}
	allocHalpha3D_ = true;

	_CalculateHalpha();
}

void Trace::CalculateHalpha(double ***halpha3d) {

	if (!_CheckHalpha()) {
		return;
	}
	
	/* allocate 1D and use pointer for 3D array */
	Halpha_ = new double[n_*nalpha_*MaxLen_];
	Halpha3D_ = halpha3d;

	allocHalpha_ = true;
	
	_CalculateHalpha();
}

void Trace::CalculateHalpha(double *halpha, double ***halpha3d) {

	if (!_CheckHalpha()) {
		return;
	}
	
	/* use pointer for both 1D and 3D arrays */
	Halpha_ = halpha;
	Halpha3D_ = halpha3d;
	
	_CalculateHalpha();
}

void Trace::_CalculateHalphaStartPoints(int i, int j,
							double *xe0, double *ye0, double *ze0,
							double *xe1, double *ye1, double *ze1) {
	
	/* this bit assumes that we are at z = 0 with vertical field lines,
	 * which will not be strictly true - this might need rewriting at some
	 * point (although I don't think it will make much of a difference */
	/* calculate the tracing start points for each alpha */
	double dt, dp, beta, dx, dy;
	
	/* dt and dp are the toroidal and poloidal components of Delta */
	dt = Delta_*cos(alpha0_[j]);
	dp = Delta_*sin(alpha0_[j]);
	
	/* rotate based on the local time */
	beta = atan2(-xfe3_[i],-yfe3_[i]);
	dx = dp*cos(beta) - dt*sin(beta);
	dy = dp*sin(beta) + dt*cos(beta);
	
	/* set the start points of the new field lines */
	xe0[0] = xfe3_[i] + dx;
	ye0[0] = yfe3_[i] + dy;
	ze0[0] = zfe3_[i];
	xe1[0] = xfe3_[i] - dx;
	ye1[0] = yfe3_[i] - dy;
	ze1[0] = zfe3_[i];

	printf("A: %f %f %f\n",xfe3_[i],yfe3_[i],zfe3_[i]);
	printf("B: %f %f %f\n",xe0[0],ye0[0],ze0[0]);
	printf("C: %f %f %f\n",xe1[0],ye1[0],ze1[0]);

}

void Trace::Field(	double x, double y, double z, 
					double *Bx, double *By, double *Bz) {
	
	/* loop through each field function, adding the field*/
	int i;
	double bx, by, bz;
	Bx[0] = 0.0;
	By[0] = 0.0;
	Bz[0] = 0.0;
	for (i=0;i<nf_;i++) {
		Funcs_[i](x,y,z,&bx,&by,&bz);
		//printf("%f %f %f \n",bx,by,bz);
		Bx[0] += bx;
		By[0] += by;
		Bz[0] += bz;
	}
}

void Trace::StepVector(	double x, double y, double z, double step3,
						double *rx, double *ry, double *rz) {
	/* based on the RHAND_08 function from GEOPACK */
	
	double bx, by, bz, s3bm;
	Field(x,y,z,&bx,&by,&bz);
	
	s3bm = step3/sqrt(bx*bx + by*by + bz*bz); 
	/* this is a unit vector scaled by 1/3 of the step size */
	rx[0] = s3bm*bx;
	ry[0] = s3bm*by;
	rz[0] = s3bm*bz;
						
}

bool Trace::ContinueTrace(double x, double y, double z, double *R) {
	
	R[0] = sqrt(x*x + y*y + z*z);
	
	/* check if we have gone too far away (this is arbitrary, I know...)*/
	if (R[0] >= MaxR_) {
		return false;
	}

	/* check that input is actually finite */
	if (!isfinite(R[0])) {
		printf("Invalid value detected\n");
		return false;
	}
	
	/* figure out some latitudes (t) */
	double rho = sqrt(x*x + y*y);
	double t = atan2(z,rho);

	/* work out the Radius of Jupiter at that latitude */
	double rhoj, zj, Rj;
	if (!SurfaceIsSphere_) {
		rhoj = as_*cos(t);
		zj = bs_*sin(t);
		Rj = sqrt(rhoj*rhoj + zj*zj);
	} else {
		Rj = rs_;
	}

	/* do the same for the ionosphere */
	double rhoi, zi, Ri;
	if (!IonosphereIsSphere_) {
		rhoi = ai_*cos(t);
		zi = bi_*sin(t);
		Ri = sqrt(rhoi*rhoi + zi*zi);
	} else {
		Ri = ri_;
	}

	/* work out the minimum value, below which we stop */	
	double rmin = std::min(Rj,Ri)*RMultiplier_;

	/* stop the trace if we go below rmin */
	if (R[0] < rmin) {
		return false;
	}
	
	return true;
	
}


void Trace::Step(	double x0, double y0, double z0,
					double *step, 
					double *x, double *y, double *z,
					double *Bx, double *By, double *Bz) {
	
	/* based on the STEP_08 function from GEOPACK */
	double rx1,ry1,rz1;
	double rx2,ry2,rz2;
	double rx3,ry3,rz3;
	double rx4,ry4,rz4;
	double rx5,ry5,rz5;
	double x1,y1,z1;
	double x2,y2,z2;
	double x3,y3,z3;
	double x4,y4,z4;
	double step3;
	double Err;
	double xn, yn, zn, Rn;
	bool cont;
	bool adjstep = false;
	bool repeat = true;	
	
	while (repeat) {
		/* this bit repeats until we get a desired step size */
		step3 = step[0]/3.0;
		StepVector(x0,y0,z0,step3,&rx1,&ry1,&rz1);
		x1 = x0 + rx1;
		y1 = y0 + ry1;
		z1 = z0 + rz1;
		StepVector(x1,y1,z1,step3,&rx2,&ry2,&rz2);
		x2 = x0 + 0.5*(rx1 + rx2);
		y2 = y0 + 0.5*(ry1 + ry2);
		z2 = z0 + 0.5*(rz1 + rz2);
		StepVector(x2,y2,z2,step3,&rx3,&ry3,&rz3);
		x3 = x0 + 0.375*(rx1 + 3*rx3);
		y3 = y0 + 0.375*(ry1 + 3*ry3);
		z3 = z0 + 0.375*(rz1 + 3*rz3);
		StepVector(x3,y3,z3,step3,&rx4,&ry4,&rz4);
		x4 = x0 + 1.5*(rx1 - 3*rx3 + 4*rx4);
		y4 = y0 + 1.5*(ry1 - 3*ry3 + 4*ry4);
		z4 = z0 + 1.5*(rz1 - 3*rz3 + 4*rz4);
		StepVector(x4,y4,z4,step3,&rx5,&ry5,&rz5);
		
		Err  = fabs(rx1 - 4.5*rx3 + 4*rx4 - 0.5*rx5);
		Err += fabs(ry1 - 4.5*ry3 + 4*ry4 - 0.5*ry5);
		Err += fabs(rz1 - 4.5*rz3 + 4*rz4 - 0.5*rz5);

		/* check that the next position is good */
		xn = x0 + 0.5*(rx1 + 4*rx4 + rx5);
		yn = y0 + 0.5*(ry1 + 4*ry4 + ry5);
		zn = z0 + 0.5*(rz1 + 4*rz4 + rz5);	
		cont = ContinueTrace(xn,yn,zn,&Rn);
		if ((!cont) && (fabs(step[0]) > MinStep_)) {
			step[0] = 0.5*step[0];
			if (fabs(step[0]) < MinStep_) {
				step[0] = (step[0]/fabs(step[0]))*MinStep_;
			}
			adjstep = true;
		} else if (!cont) {
			step[0] = (step[0]/fabs(step[0]))*MinStep_;
			repeat = false;
		}
		
		if (cont) {
			if ((Err <= ErrMax_) && (fabs(step[0]) <= MaxStep_)) {
				repeat = false;
			} else {
				if (Err > ErrMax_) {
					if (fabs(step[0]) > MinStep_) {
						step[0] = step[0]*0.5;
					} else {
						repeat = false;
					}
				}
				if (fabs(step[0]) > MaxStep_) {
					step[0] = (step[0]/fabs(step[0]))*MaxStep_;
				}
			}
			
			if ((Err < 0.04*ErrMax_) && (fabs(step[0]) < (MaxStep_/1.5))) {
				step[0] = 1.5*step[0];
			}		
		}
	}
	
	x[0] = xn;
	y[0] = yn;
	z[0] = zn;	
	
	Field(x[0],y[0],z[0],Bx,By,Bz);					
}

void Trace::ReverseElements(int n, double *x) {
	int i;
	double tmp;
	for (i=0;i<(n/2);i++) {
		tmp = x[i];
		x[i] = x[n-i-1];
		x[n-i-1] = tmp;
	}
}

void Trace::RKMTrace(	double x0, double y0, double z0,
						int *nstep, double *R,
						double *x, double *y, double *z,
						double *Bx, double *By, double *Bz) {

	/* intialize the trace */
	nstep[0] = 1;
	x[0] = x0;
	y[0] = y0;
	z[0] = z0;
	Field(x0,y0,z0,&Bx[0],&By[0],&Bz[0]);

	double step;
	bool cont = ContinueTrace(x[0],y[0],z[0],&R[0]);
	
	/* trace in one direction */
	if ((TraceDir_ == 1) || (TraceDir_ == 0)) {
		/* I think this will trace opposite to the direction of the field,
		 * into the northern hemisphere */ 
		step = -InitStep_;
		while ((cont) && (nstep[0] < (MaxLen_/2 - 1))) {
			Step(	x[nstep[0]-1],y[nstep[0]-1],z[nstep[0]-1],&step,
					&x[nstep[0]],&y[nstep[0]],&z[nstep[0]],
					&Bx[nstep[0]],&By[nstep[0]],&Bz[nstep[0]]);
			cont = ContinueTrace(x[nstep[0]],y[nstep[0]],z[nstep[0]],&R[nstep[0]]);
			nstep[0]++;
		}
	}
	
	/* reverse the elements of the trace */
	ReverseElements(nstep[0],x);
	ReverseElements(nstep[0],y);
	ReverseElements(nstep[0],z);
	ReverseElements(nstep[0],Bx);
	ReverseElements(nstep[0],By);
	ReverseElements(nstep[0],Bz);
	ReverseElements(nstep[0],R);
	
	/* trace in the opposite direction */
	cont = ContinueTrace(x[nstep[0]-1],y[nstep[0]-1],z[nstep[0]-1],&R[nstep[0]-1]);
	if ((TraceDir_ == -1) || (TraceDir_ == 0)) {
		/* hopefully this will go in the direction fo the field vectors
		 * towards the southern hemisphere */
		step = InitStep_;
		while ((cont) && (nstep[0] < (MaxLen_ - 1))) {
			Step(	x[nstep[0]-1],y[nstep[0]-1],z[nstep[0]-1],&step,
					&x[nstep[0]],&y[nstep[0]],&z[nstep[0]],
					&Bx[nstep[0]],&By[nstep[0]],&Bz[nstep[0]]);
			cont = ContinueTrace(x[nstep[0]],y[nstep[0]],z[nstep[0]],&R[nstep[0]]);
			nstep[0]++;
		}
	}
	
	/* sort the footprints out */
	//FixFootprints(nstep[0],R,x,y,z,Bx,By,Bz);
}

void Trace::FixFootprints(	int nstep, double *R,
						double *x, double *y, double *z,
						double *Bx, double *By, double *Bz) {
	
	/* this function will fix each of the footprints using interpolation
	 * such that they were actually on the surface */
	
	/* these are temporary variables to be used for interpolation */
	int i0, i1;
	double a = 1.0;
	double b = 0.935;	
	double rho[2], t[2], rhoj[2], zj[2];
	double ts[2], tr[2], Rs, rs[2];
	double dx, dy, dz, dr, ds, s1, m;
	
	if ((TraceDir_ == 1) || (TraceDir_ == 0)) {
		/* north footprint */
		i0 = 0;
		i1 = 1;
		
		/* get the approximate R where the "surface" is */
		rho[0] = sqrt(x[i0]*x[i0] + y[i0]*y[i0]);
		rho[1] = sqrt(x[i1]*x[i1] + y[i1]*y[i1]);
		t[0] = atan2(z[i0],rho[0]);
		t[1] = atan2(z[i1],rho[1]);
		rhoj[0] = a*cos(t[0]);
		rhoj[1] = a*cos(t[1]);
		zj[0] = b*sin(t[0]);
		zj[1] = b*sin(t[1]);
		rs[0] = sqrt(rhoj[0]*rhoj[0] + zj[0]*zj[0]);
		rs[1] = sqrt(rhoj[1]*rhoj[1] + zj[1]*zj[1]);
		Rs = 0.5*(rs[0] + rs[1]);
		
		/* use linear interpolation (it's generally close enough) */
		dx = x[i1] - x[i0];
		dy = y[i1] - y[i0];
		dz = z[i1] - z[i0];
		ts[0] = 0.0;
		ts[1] = sqrt(dx*dx + dy*dy + dz*dz);
		tr[0] = sqrt(x[i0]*x[i0] + y[i0]*y[i0] + z[i0]*z[i0]);
		tr[1] = sqrt(x[i1]*x[i1] + y[i1]*y[i1] + z[i1]*z[i1]);
			
		/* gradient */
		dr = tr[1] - tr[0];
		ds = ts[1] - ts[0];
		m = dr/ds;
			
		/* solve for R = 1, where R = m*s + c */
		s1 = (Rs - tr[0])/m;
			
		/* interpolate x, y and z*/
		x[i0] = (dx/ds)*s1 + x[i0];
		y[i0] = (dy/ds)*s1 + y[i0];
		z[i0] = (dz/ds)*s1 + z[i0];
		
		/* recalculate the field */
		Field(x[i0],y[i0],z[i0],&Bx[i0],&By[i0],&Bz[i0]);	
		
		/* fix R */
		R[i0] = sqrt(x[i0]*x[i0] + y[i0]*y[i0] + z[i0]*z[i0]);
	}
	
	if ((TraceDir_ == -1) || (TraceDir_ == 0)) {
		/* south footprint */
		i0 = nstep - 1;
		i1 = nstep - 2;

		/* get the approximate R where the "surface" is */
		rho[0] = sqrt(x[i0]*x[i0] + y[i0]*y[i0]);
		rho[1] = sqrt(x[i1]*x[i1] + y[i1]*y[i1]);
		t[0] = atan2(z[i0],rho[0]);
		t[1] = atan2(z[i1],rho[1]);
		rhoj[0] = a*cos(t[0]);
		rhoj[1] = a*cos(t[1]);
		zj[0] = b*sin(t[0]);
		zj[1] = b*sin(t[1]);
		rs[0] = sqrt(rhoj[0]*rhoj[0] + zj[0]*zj[0]);
		rs[1] = sqrt(rhoj[1]*rhoj[1] + zj[1]*zj[1]);
		Rs = 0.5*(rs[0] + rs[1]);
		
		/* use linear interpolation (it's generally close enough) */
		dx = x[i1] - x[i0];
		dy = y[i1] - y[i0];
		dz = z[i1] - z[i0];
		ts[0] = 0.0;
		ts[1] = sqrt(dx*dx + dy*dy + dz*dz);
		tr[0] = sqrt(x[i0]*x[i0] + y[i0]*y[i0] + z[i0]*z[i0]);
		tr[1] = sqrt(x[i1]*x[i1] + y[i1]*y[i1] + z[i1]*z[i1]);
			
		/* gradient */
		dr = tr[1] - tr[0];
		ds = ts[1] - ts[0];
		m = dr/ds;
			
		/* solve for R = 1, where R = m*s + c */
		s1 = (Rs - tr[0])/m;
			
		/* interpolate x, y and z*/
		x[i0] = (dx/ds)*s1 + x[i0];
		y[i0] = (dy/ds)*s1 + y[i0];
		z[i0] = (dz/ds)*s1 + z[i0];

		/* recalculate the field */
		Field(x[i0],y[i0],z[i0],&Bx[i0],&By[i0],&Bz[i0]);	
		
		/* fix R */
		R[i0] = sqrt(x[i0]*x[i0] + y[i0]*y[i0] + z[i0]*z[i0]);
	}
}


void Trace::TraceField(	int *nstep,
						double **x, double **y, double **z, double **R,
						double **bx, double **by, double **bz) {
	
	/* link the pointers within the object to those supplied by this 
	 * function					*/
	nstep_ = nstep;
	x_ = x;					
	y_ = y;					
	z_ = z;					
	bx_ = bx;					
	by_ = by;					
	bz_ = bz;	
	R_ = R;	

	/* call the tracing code */
	_TraceField();
}


void Trace::TraceField() {
	
	/* no pointers provided: allocate them*/
	nstep_ = new int[n_];
	x_ = new double*[n_];					
	y_ = new double*[n_];					
	z_ = new double*[n_];					
	bx_ = new double*[n_];					
	by_ = new double*[n_];					
	bz_ = new double*[n_];
	R_ = new double*[n_];
	int i;
	for (i=0;i<n_;i++) {
		x_[i] = new double[MaxLen_];					
		y_[i] = new double[MaxLen_];					
		z_[i] = new double[MaxLen_];					
		bx_[i] = new double[MaxLen_];					
		by_[i] = new double[MaxLen_];					
		bz_[i] = new double[MaxLen_];		
		R_[i] = new double[MaxLen_];		
	}		
	allocTrace_ = true;
	
	/* call the tracing code */
	_TraceField();
	
}


void Trace::_TraceField() {
	
	/* this function actually calls the tracing routines */
	/* check this hasn't already been done */
	if (tracedField_) {
		printf("Attempted to trace twice? not happening mate...\n");
		return;
	}
	
	/* check we have input positions */
	if (!inputPos_) {
		printf("Need InputPos() before trace\n");
		return;
	}

	int i;
	for (i=0;i<n_;i++) {
		if (Verbose_) {
			printf("\rTracing field line %d of %d (%6.2f)%%",i+1,n_,((float) (i+1)*100.0)/n_);
		}


		/* perform trace */
		RKMTrace(x0_[i],y0_[i],z0_[i],&nstep_[i],R_[i],x_[i],y_[i],z_[i],bx_[i],by_[i],bz_[i]);
		
		
	}	
	if (Verbose_) { 
		printf("\n");
	}
	tracedField_ = true;
}


void Trace::CalculateTraceDist() {
	int i;
	S_ = new double*[n_];
	for (i=0;i<n_;i++) {
		S_[i] = new double[MaxLen_];
	}
	allocDist_ = true;
	
	_CalculateTraceDist();
}

void Trace::CalculateTraceDist(double **S) {

	S_ = S;
	_CalculateTraceDist();
}


void Trace::_CalculateTraceDist() {
	int i, j;
	double dx, dy, dz;
	for (i=0;i<n_;i++) {
		S_[i][0] = 0.0;
		for (j=1;j<nstep_[i];j++) {
			dx = x_[i][j] - x_[i][j-1];
			dy = y_[i][j] - y_[i][j-1];
			dz = z_[i][j] - z_[i][j-1];
			S_[i][j] = S_[i][j-1] + sqrt(dx*dx + dy*dy + dz*dz);
		}
	}
	hasDist_ = true;
}



void Trace::CalculateTraceRnorm() {
	int i;
	Rnorm_ = new double*[n_];
	for (i=0;i<n_;i++) {
		Rnorm_[i] = new double[MaxLen_];
	}
	allocRnorm_ = true;
	
	_CalculateTraceRnorm();
}

void Trace::CalculateTraceRnorm(double **Rnorm) {
	
	Rnorm_ = Rnorm;
	
	_CalculateTraceRnorm();
}

void Trace::_CalculateTraceRnorm() {
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<nstep_[i];j++) {
			/* need footprints and R done first */
			
		}
	}
	hasRnorm_ = true;
}


void Trace::CalculateTraceFP() {
	int i;
	FP_ = new double*[n_];
	for (i=0;i<n_;i++) {
		FP_[i] = new double[20];
	}
	allocFootprints_ = true;
	
	_CalculateTraceFP();
}

void Trace::CalculateTraceFP(double **FP) {
	
	FP_ = FP;
	
	_CalculateTraceFP();
}

void Trace::_CalculateTraceFP() {

	/* before running this we should check that the following functions 
	 * have been called:
	 * 1. TraceField()
	 * 2. CalcualteTraceDist() */
	if (!tracedField_) {
		printf("Call TraceField() before calculating footprints\n");
		return;
	}
	if (!hasDist_) {
		printf("Call CalcualteTraceDist() before calculating footprints\n");
		return;
	}
	

	/* allocate the surface footprints (SIII) */
	xfn3_ = new double[n_];
	yfn3_ = new double[n_];
	zfn3_ = new double[n_];
	xfs3_ = new double[n_];
	yfs3_ = new double[n_];
	zfs3_ = new double[n_];

	/* allocate ionospheric footprints (SIII) */
	xin3_ = new double[n_];
	yin3_ = new double[n_];
	zin3_ = new double[n_];
	xis3_ = new double[n_];
	yis3_ = new double[n_];
	zis3_ = new double[n_];
	
	/* and "equatorial" footprints" (SIII) */	
	xfe3_ = new double[n_];
	yfe3_ = new double[n_];
	zfe3_ = new double[n_];

	
	/* allocate the surface footprints (Mag) */
	xfnm_ = new double[n_];
	yfnm_ = new double[n_];
	zfnm_ = new double[n_];
	xfsm_ = new double[n_];
	yfsm_ = new double[n_];
	zfsm_ = new double[n_];

	/* allocate ionospheric footprints (Mag) */
	xinm_ = new double[n_];
	yinm_ = new double[n_];
	zinm_ = new double[n_];
	xism_ = new double[n_];
	yism_ = new double[n_];
	zism_ = new double[n_];
	
	/* and "equatorial" footprints" (Mag) */	
	xfem_ = new double[n_];
	yfem_ = new double[n_];
	zfem_ = new double[n_];
	
	allocEqFP_ = true;
	
	double rho, latn, lats, lonn, lons, lone, Lshell, FlLen;
	double mlatn, mlats, mlonn, mlons, mlone;
	double latni, latsi, lonni, lonsi;
	double mlatni, mlatsi, mlonni, mlonsi;

	int i, j, imaxR;
	for (i=0;i<n_;i++) {
		printf("Footprint %d of %d\n",i+1,n_);

		/* initialize with NAN */
		xfn3_[i] = NAN;
		yfn3_[i] = NAN;
		zfn3_[i] = NAN;

		xfs3_[i] = NAN;
		yfs3_[i] = NAN;
		zfs3_[i] = NAN;

		xin3_[i] = NAN;
		yin3_[i] = NAN;
		zin3_[i] = NAN;

		xis3_[i] = NAN;
		yis3_[i] = NAN;
		zis3_[i] = NAN;

		xfnm_[i] = NAN;
		yfnm_[i] = NAN;
		zfnm_[i] = NAN;

		xfsm_[i] = NAN;
		yfsm_[i] = NAN;
		zfsm_[i] = NAN;

		xinm_[i] = NAN;
		yinm_[i] = NAN;
		zinm_[i] = NAN;

		xism_[i] = NAN;
		yism_[i] = NAN;
		zism_[i] = NAN;

		xfe3_[i] = NAN;
		yfe3_[i] = NAN;
		zfe3_[i] = NAN;

		xfem_[i] = NAN;
		yfem_[i] = NAN;
		zfem_[i] = NAN;



		/* calculate the surface footprints */
		if (SurfaceIsSphere_) {
			printf("Sphere Surface\n");
			footprints(nstep_[i],x_[i],y_[i],z_[i],rs_,rs_,xt_,xp_,
						&xfn3_[i],&yfn3_[i],&zfn3_[i],&xfs3_[i],&yfs3_[i],&zfs3_[i]);
		} else {
			printf("Spheroid surface\n");
			footprints(nstep_[i],x_[i],y_[i],z_[i],as_,bs_,xt_,xp_,
						&xfn3_[i],&yfn3_[i],&zfn3_[i],&xfs3_[i],&yfs3_[i],&zfs3_[i]);
			printf("S %f %f %f %f %f %f\n",xfn3_[i],yfn3_[i],zfn3_[i],xfs3_[i],yfs3_[i],zfs3_[i]);
		}

		/* calculate the surface footprints */
		if (IonosphereIsSphere_) {
			printf("Sphere Ionosphere\n");
			footprints(nstep_[i],x_[i],y_[i],z_[i],ri_,ri_,xt_,xp_,
						&xin3_[i],&yin3_[i],&zin3_[i],&xis3_[i],&yis3_[i],&zis3_[i]);
			printf("I %f\n",xin3_[i]);
			printf("I %f\n",yin3_[i]);
			printf("I %f\n",zin3_[i]);
			printf("I %f\n",xis3_[i]);
			printf("I %f\n",yis3_[i]);
			printf("I %f\n",zis3_[i]);

		} else {
			printf("Spheroid Ionosphere\n");
			footprints(nstep_[i],x_[i],y_[i],z_[i],ai_,bi_,xt_,xp_,
						&xin3_[i],&yin3_[i],&zin3_[i],&xis3_[i],&yis3_[i],&zis3_[i]);
		}

		if (!isnan(xfn3_[i]) && !isnan(xfs3_[i])) {
			printf("Equator\n");
			eqfootprints(nstep_[i],x_[i],y_[i],z_[i],&xfe3_[i],&yfe3_[i],&zfe3_[i],&Lshell,&lone);
			FlLen = S_[i][nstep_[i]];
		} else {
			Lshell = NAN;
			lone = NAN;
			FlLen = NAN;
		}


		/* convert to magnetic coordinates */
		if (!isnan(xfn3_[i])) {
			printf("Mag Conversion North\n");
			SIIItoMag(xfn3_[i],yfn3_[i],zfn3_[i],xt_,xp_,&xfnm_[i],&yfnm_[i],&zfnm_[i]);

		 	rho = sqrt(xfn3_[i]*xfn3_[i] + yfn3_[i]*yfn3_[i]);
			
		 	latn = rad2deg*atan2(zfn3_[i],rho);
		 	lonn = rad2deg*atan2(yfn3_[i],xfn3_[i]);
		 	printf("%f %f %f\n",zfn3_[i],rho,latn);
			printf("%f %f %f\n",yfn3_[i],xfn3_[i],lonn);
			rho = sqrt(xfnm_[i]*xfnm_[i] + yfnm_[i]*yfnm_[i]);
		 	mlatn = rad2deg*atan2(zfnm_[i],rho);
		 	mlonn = rad2deg*atan2(yfnm_[i],xfnm_[i]);
			printf("%f %f %f\n",zfnm_[i],rho,mlatn);
			printf("%f %f %f\n",yfnm_[i],xfnm_[i],mlonn);
		} else {
			latn = NAN;
			lonn = NAN;

			mlatn = NAN;
			mlonn = NAN;

		}
		if (!isnan(xin3_[i])) {
			SIIItoMag(xin3_[i],yin3_[i],zin3_[i],xt_,xp_,&xinm_[i],&yinm_[i],&zinm_[i]);

		 	rho = sqrt(xin3_[i]*xin3_[i] + yin3_[i]*yin3_[i]);
		 	latni = rad2deg*atan2(zin3_[i],rho);
		 	lonni = rad2deg*atan2(yin3_[i],xin3_[i]);
		 	printf("%f %f %f\n",zin3_[i],rho,latni);
			printf("%f %f %f\n",yin3_[i],xin3_[i],lonni);
		 	
			rho = sqrt(xinm_[i]*xinm_[i] + yinm_[i]*yinm_[i]);
		 	mlatni = rad2deg*atan2(zinm_[i],rho);
		 	mlonni = rad2deg*atan2(yinm_[i],xinm_[i]);
			printf("%f %f %f\n",zinm_[i],rho,mlatni);
			printf("%f %f %f\n",yinm_[i],xinm_[i],mlonni);
		} else {
			latni = NAN;
			lonni = NAN;

			mlatni = NAN;
			mlonni = NAN;
		}
		if (!isnan(xfs3_[i])) {
			SIIItoMag(xfs3_[i],yfs3_[i],zfs3_[i],xt_,xp_,&xfsm_[i],&yfsm_[i],&zfsm_[i]);

		 	rho = sqrt(xfs3_[i]*xfs3_[i] + yfs3_[i]*yfs3_[i]);
		 	lats = rad2deg*atan2(zfs3_[i],rho);
		 	lons = rad2deg*atan2(yfs3_[i],xfs3_[i]);
		 	
			rho = sqrt(xfsm_[i]*xfsm_[i] + yfsm_[i]*yfsm_[i]);
		 	mlats = rad2deg*atan2(zfsm_[i],rho);
		 	mlons = rad2deg*atan2(yfsm_[i],xfsm_[i]);

		} else {
			lats = NAN;
			lons = NAN;

			mlats = NAN;
			mlons = NAN;
		}
		if (!isnan(xis3_[i])) {
			SIIItoMag(xis3_[i],yis3_[i],zis3_[i],xt_,xp_,&xism_[i],&yism_[i],&zism_[i]);

		 	rho = sqrt(xis3_[i]*xis3_[i] + yis3_[i]*yis3_[i]);
		 	latsi = rad2deg*atan2(zis3_[i],rho);
		 	lonsi = rad2deg*atan2(yis3_[i],xis3_[i]);
		 	
			rho = sqrt(xism_[i]*xism_[i] + yism_[i]*yism_[i]);
		 	mlatsi = rad2deg*atan2(zism_[i],rho);
		 	mlonsi = rad2deg*atan2(yism_[i],xism_[i]);
		} else {

			latsi = NAN;
			lonsi = NAN;

			mlatsi = NAN;
			mlonsi = NAN;
		}
		if (!isnan(xfe3_[i])) {
			SIIItoMag(xfe3_[i],yfe3_[i],zfe3_[i],xt_,xp_,&xfem_[i],&yfem_[i],&zfem_[i]);
			mlone = rad2deg*atan2(yfem_[i],xfem_[i]);
		}


		/* calculate the equatorial footprints */

		// /* north footprint */
		// if ((TraceDir_ == 0) || (TraceDir_ == 1)) {
		// 	xfn_[i] = x_[i][0];
		// 	yfn_[i] = y_[i][0];
		// 	zfn_[i] = z_[i][0];
			
		// 	/* latitude (not co-lat)  and longitude - in degrees*/
		// 	rho = sqrt(xfn_[i]*xfn_[i] + yfn_[i]*yfn_[i]);
		// 	latn = rad2deg*atan2(zfn_[i],rho);
		// 	lonn = rad2deg*atan2(yfn_[i],xfn_[i]);
		// } else {
		// 	latn = NAN;
		// 	lonn = NAN;
		// }
		
		// /* south footprint */
		// if ((TraceDir_ == 0) || (TraceDir_ == -1)) {
		// 	xfs_[i] = x_[i][nstep_[i]-1];
		// 	yfs_[i] = y_[i][nstep_[i]-1];
		// 	zfs_[i] = z_[i][nstep_[i]-1];
			
		// 	/* latitude (not co-lat)  and longitude - in degrees*/
		// 	rho = sqrt(xfs_[i]*xfs_[i] + yfs_[i]*yfs_[i]);
		// 	lats = rad2deg*atan2(zfs_[i],rho);
		// 	lons = rad2deg*atan2(yfs_[i],xfs_[i]);
		// } else {
		// 	lats = NAN;
		// 	lons = NAN;
		// }				
		
		// /* equatorial (sort of, axtually at Rmax) footprint */
		// if (TraceDir_ == 0) {
		// 	/* find the furthest point along the field line */
		// 	imaxR = -1;
		// 	Lshell = 0.0;
		// 	for (j=0;j<nstep_[i];j++) {
		// 		if (R_[i][j] > Lshell) {
		// 			Lshell = R_[i][j];
		// 			imaxR = j;
		// 		}
		// 	}
			
		// 	xfe_[i] = x_[i][imaxR];
		// 	yfe_[i] = y_[i][imaxR];
		// 	zfe_[i] = z_[i][imaxR];
			
		// 	/* latitude (not co-lat)  and longitude - in degrees*/
		// 	lone = rad2deg*atan2(yfe_[i],xfe_[i]);
			
		// 	/* field length */
		// 	FlLen = S_[i][nstep_[i]];
		// } else {
		// 	Lshell = NAN;
		// 	lone = NAN;
		// 	FlLen = NAN;
		// }	
		FP_[i][0] = latn;
		FP_[i][1] = lonn;
		FP_[i][2] = lats;
		FP_[i][3] = lons;
		FP_[i][4] = latni;
		FP_[i][5] = lonni;
		FP_[i][6] = latsi;
		FP_[i][7] = lonsi;
		FP_[i][8] = lone;
		FP_[i][9] = mlatn;
		FP_[i][10] = mlonn;
		FP_[i][11] = mlats;
		FP_[i][12] = mlons;
		FP_[i][13] = mlatni;
		FP_[i][14] = mlonni;
		FP_[i][15] = mlatsi;
		FP_[i][16] = mlonsi;
		FP_[i][17] = mlone;
		FP_[i][18] = Lshell;
		FP_[i][19] = FlLen;
	}
	hasFootprints_ = true;
}



void Trace::GetTrace(double **x,double **y, double **z) {
	/* copy position into output arrays*/
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<nstep_[i];j++) {
			x[i][j] = x_[i][j];
			y[i][j] = y_[i][j];
			z[i][j] = z_[i][j];
		}
	}
}

void Trace::GetTrace(	double **x,double **y, double **z,
						double **Bx,double **By, double **Bz) {
	/* copy field into output arrays*/
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<nstep_[i];j++) {
			Bx[i][j] = bx_[i][j];
			By[i][j] = by_[i][j];
			Bz[i][j] = bz_[i][j];
		}
	}
	
	/* get the position */
	GetTrace(x,y,z);
}


void Trace::GetTraceDist(double **S) {
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<nstep_[i];j++) {
			S[i][j] = S_[i][j];
		}
	}	
}

void Trace::GetTraceR(double **R) {
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<nstep_[i];j++) {
			R[i][j] = R_[i][j];
		}
	}	
}

void Trace::GetTraceRnorm(double **Rnorm) {
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<nstep_[i];j++) {
			Rnorm[i][j] = Rnorm_[i][j];
		}
	}	
}

void Trace::GetTraceFootprints(double **FP) {
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<7;j++) {
			FP[i][j] = FP_[i][j];
		}
	}	
}

void Trace::GetTraceNstep(int *nstep) {
	int i;
	for (i=0;i<n_;i++) {
		nstep[i] = nstep_[i];
	}
}


void Trace::SetTraceBoundDefaults() {

	/* set the default radius to 1.0 */
	rs_ = 1.0;
	as_ = 1.0;
	bs_ = 0.93513;
	SurfaceIsSphere_ = false;

	/* ionoshere = surface + 500 km,
	although here we use the polar radius (66854 km)
	so ri_ = (66854 + 500)/71492*/
	ri_ = 0.94212;
	ai_ = 1.00699;
	bi_ = 0.94212;
	IonosphereIsSphere_ = true;

	/* maximum radial distance */
	MaxR_ = 1000.0;

}

void Trace::SetSurfaceSpheroidR(double a, double b) {

	as_ = a;
	bs_ = b;
}

void Trace::GetSurfaceSpheroidR(double *a, double *b) {

	*a = as_;
	*b = bs_;
}

void Trace::SetSurfaceSphereR(double r) {

	rs_ = r;
}

double Trace::GetSurfaceSphereR() {

	return rs_;
}

void Trace::SetSurfaceIsSphere(bool s) {

	SurfaceIsSphere_ = s;
}

bool Trace::GetSurfaceIsSphere() {

	return SurfaceIsSphere_;
}


void Trace::SetIonosphereSpheroidR(double a, double b) {

	ai_ = a;
	bi_ = b;
}

void Trace::GetIonosphereSpheroidR(double *a, double *b) {

	*a = ai_;
	*b = bi_;
}

void Trace::SetIonosphereSphereR(double r) {

	ri_ = r;
}

double Trace::GetIonosphereSphereR() {

	return ri_;
}

void Trace::SetIonosphereIsSphere(bool s) {

	IonosphereIsSphere_ = s;
}

bool Trace::GetIonosphereIsSphere() {

	return IonosphereIsSphere_;
}

void Trace::SetTraceMaxR(double MaxR) {

	MaxR_ = MaxR;
}

double Trace::GetTraceMaxR() {

	return MaxR_;
}




void Trace::SetMagTilt(double xt) {
	/* input in degrees */
	xt_ = xt*deg2rad;
}

double Trace::GetMagTilt() {
	/* convert back to degrees */
	return xt_*rad2deg;
}

void Trace::SetMagTiltAzimuth(double xp) {
	xp_ = xp*deg2rad;
}

double Trace::GetMagTiltAzimuth() {
	return xp_*rad2deg;
}
	