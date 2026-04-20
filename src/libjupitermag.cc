#include "libjupitermag.h"

extern "C" void JupitermagGetCon2020Params(
	double *mui, double *irho, double *r0, double *r1,
	double *d, double *xt, double *xp, char *eqtype,
	bool *Edwards, bool *ErrChk, bool *CartIn, bool *CartOut,
	bool *smooth, double *DeltaRho, double *DeltaZ,
	double *g, char *azfunc, double *wO_open, double *wO_om,
	double *thetamm, double *dthetamm, double *thetaoc, double *dthetaoc) {
	con2020::GetCon2020Params(mui, irho, r0, r1, d, xt, xp, eqtype,
					 Edwards, ErrChk, CartIn, CartOut,
					 smooth, DeltaRho, DeltaZ, g, azfunc,
					 wO_open, wO_om, thetamm, dthetamm, thetaoc, dthetaoc);
}

extern "C" void JupitermagSetCon2020Params(
	double mui, double irho, double r0, double r1,
	double d, double xt, double xp, const char *eqtype,
	bool Edwards, bool ErrChk, bool CartIn, bool CartOut,
	bool smooth, double DeltaRho, double DeltaZ,
	double g, const char *azfunc, double wO_open, double wO_om,
	double thetamm, double dthetamm, double thetaoc, double dthetaoc) {
	con2020::SetCon2020Params(mui, irho, r0, r1, d, xt, xp, eqtype,
					 Edwards, ErrChk, CartIn, CartOut,
					 smooth, DeltaRho, DeltaZ, g, azfunc,
					 wO_open, wO_om, thetamm, dthetamm, thetaoc, dthetaoc);
}

extern "C" void JupitermagSetInternalCFG(const char *Model, bool CartIn, bool CartOut, int MaxDeg) {
	SetInternalCFG(Model, CartIn, CartOut, MaxDeg);
}

extern "C" void JupitermagGetInternalCFG(char *Model, bool *CartIn, bool *CartOut, int *MaxDeg) {
	GetInternalCFG(Model, CartIn, CartOut, MaxDeg);
}


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
				int nalpha, double *alpha, double *halpha) {
	
	/* before calling this wrapper function, any field models used should
	 * be configured. This function will not do any of that, so strange
	 * things could happen. Make sure that all models are Cartesian in 
	 * and out! */
	std::vector<jupitermag::FieldFuncPtr> Funcs;

	/* internal model */
	Funcs.push_back(getModelFieldPtr(IntFunc));

	/* external model */
	int i;
	for (i=0;i<nExt;i++) {
		if (strcmp(ExtFunc[i],"Con2020") == 0) {
			Funcs.push_back(con2020::Con2020Field);
		}
	}

	/* if there are no functions then return */
	if (Funcs.size() == 0) {
		printf("No valid model functions provided\n");
		return false;
	}

	/* initialise the trace object */
	jupitermag::Trace T(Funcs);

	/* add the starting positions fo the traces */
	T.InputPos(n,x0,y0,z0);

	/* configure the trace parameters */
	T.SetTraceCFG(MaxLen,MaxStep,InitStep,MinStep,ErrMax,Delta,Verbose,TraceDir);

	/* set up the surfaces */
	if (ai == bi) {
		T.SetIonosphereIsSphere(true);
		T.SetIonosphereSphereR(ai);
	} else {
		T.SetIonosphereIsSphere(false);
		T.SetIonosphereSpheroidR(ai,bi);
	}
	if (as == bs) {
		T.SetSurfaceIsSphere(true);
		T.SetSurfaceSphereR(as);
	} else {
		T.SetSurfaceIsSphere(false);
		T.SetSurfaceSpheroidR(as,bs);
	}

	/* set up the alpha calculation */
	if (nalpha > 0) {
		T.SetAlpha(nalpha,alpha);
	}

	/* Trace */
	T.TraceField(nstep,x,y,z,R,Bx,By,Bz,traceRegion);

	/* trace distance, footprints, Rnorm */
	if (TraceDir == 0) {
		T.CalculateTraceDist(S);
		T.CalculateTraceFP(FP);
		T.CalculateTraceRnorm(Rnorm);
	}

	/* halpha */
	if ((nalpha > 0) && (TraceDir == 0)) {
		T.CalculateHalpha(halpha);
	}

	return true;
}
