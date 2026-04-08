
#ifndef __LIBJUPITERMAG_H__
#define __LIBJUPITERMAG_H__
#define _USE_MATH_DEFINES

#ifdef __cplusplus
	#include <algorithm>
	#include <map>
	#include <math.h>
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include <string>
	#include <tuple>
	#include <vector>
#else
	#include <math.h>
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include <stdbool.h>
#endif

#include "con2020.h"
#include "spline.h"
#include "internalfield.h"


#define LIBJUPITERMAG_VERSION_MAJOR 1
#define LIBJUPITERMAG_VERSION_MINOR 5
#define LIBJUPITERMAG_VERSION_PATCH 0


#ifdef __cplusplus
extern "C" {
#endif

	void ModelField(double p0, double p1, double p2, 
					const char *internal, const char *external, 
					bool CartIn, bool CartOut,
					double *B0, double *B1, double *B2);

	void ModelFieldArray(	int n, double *p0, double *p1, double *p2, 
							const char *internal, const char *external, 
							bool CartIn, bool CartOut,
							double *B0, double *B1, double *B2);
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
#ifdef __cplusplus
} /* extern "C" */

namespace jupitermag {


void MagtoSIII(	double xm, double ym, double zm, 
				double xt, double xp,
				double *x3, double *y3, double *z3);

void SIIItoMag(	double x3, double y3, double z3, 
				double xt, double xp,
				double *xm, double *ym, double *zm);


typedef struct FPstr {
	/* SIII North */
	double xn3;
	double yn3;
	double zn3;
	/* SIII South */
	double xs3;
	double ys3;
	double zs3;
	/* MAG North */
	double xnm;
	double ynm;
	double znm;
	/* MAG South */
	double xsm;
	double ysm;
	double zsm;

	/* SIII North */
	double lonn;
	double latn;
	/* MAG North */
	double mlonn;
	double mlatn;
	/* SIII South */
	double lons;
	double lats;
	/* MAG South */
	double mlons;
	double mlats;

} FPstr;

typedef struct EqFPstr {
	/* SIII coords */
	double x3;
	double y3;
	double z3;

	/* MAG coords */
	double xm;
	double ym;
	double zm;

	/*Equatorial footprint*/
	double lshell;
	double mlone;
	double fllen;
} EqFPstr;


void FillFPOutputArray(int n, FPstr *fpi, FPstr *fps, EqFPstr *fpe, double **FP);

void _fillfp(FPstr *fp, double xt, double xp);


void _getbegfp(	int n, double *x, double *y, double *z,
				bool begfp, bool endfp, int indmxr,
				double a, double b,
				double *xfp, double *yfp, double *zfp);

void _getendfp(	int n, double *x, double *y, double *z,
				bool begfp, bool endfp, int indmxr,
				double a, double b,
				double *xfp, double *yfp, double *zfp);

bool _posisfp(double x, double y, double z, double a, double b);

int _maxR(double n, double *x, double *y, double *z);

void _nsends(int n, double *x, double *y, double *z,
			double xt, double xp,
			bool begfp, bool endfp,
			int *begns, int *endns, int *indmxr);

double _fllen(int n, double *x, double *y, double *z) ;

void calculateEquatorialFootprints(int n, double *x, double *y, double *z,
						double xt, double xp, EqFPstr *efp);

void calculateFootprints(int n, double *x, double *y, double *z,
							double a, double b, double xt, double xp,
							FPstr *fp);

void footprints(int n, double *x, double *y, double *z,
				double a, double b, double xt, double xp,
				double *xfn, double *yfn, double *zfn,
				double *xfs, double *yfs, double *zfs);

double planetRadius(double x, double y, double z, double a, double b);

double thetaPlanetRadius(double theta, double a, double b);

bool isCrossing(double x0, double y0, double z0,
				double x1, double y1, double z1,
				double a, double b);

void interpCrossing(double x0, double y0, double z0,
					double x1, double y1, double z1,
					double a, double b,
					double *xfp, double *yfp, double *zfp);

void findFootprint(	double *x, double *y, double *z,
					int starti, int endi,
					double a, double b,
					double *xfp, double *yfp, double *zfp);

void eqfootprints(	int n, double *x, double *y, double *z,
					double *xfe, double *yfe, double *zfe,
					double *L, double *Lon);




void interptraceClosestPos(	int n, double *x, double *y, double *z,
						double *bx, double *by, double *bz,
						int n0, double *x0, double *y0, double *z0, double *s0,
						int n1, double *x1, double *y1, double *z1, double *s1,
						double *xc0, double *yc0, double *zc0,
						double *xc1, double *yc1, double *zc1 );
						
double ClosestS(double x, double y, double z,
				int nt, double *xt, double *yt, double *zt,
				double *st);
				
double AngleDiff( 	double s,								/* current position along the field line */
					spline::Spline Sx, spline::Spline Sy, spline::Spline Sz,	/* Splines converting s to a  vector */
					double x, double y, double z,		/* this is the position along the original field line */
					double bx, double by, double bz);
					
					
void OptimizePos(	double x, double y, double z,
					double bx, double by, double bz,
					double s0, 
					spline::Spline Sx, spline::Spline Sy, spline::Spline Sz,
					double *xc, double *yc, double *zc);


typedef void (*FieldFuncPtr)(double,double,double,double*,double*,double*);

typedef std::tuple<bool,int> BoolIntTuple;

class Trace {
	
	public:
		Trace(std::vector<FieldFuncPtr>);
		~Trace();
		
		void InputPos(int,double*,double*,double*);
		void SetTraceCFG(int, double,double,double,double,double,bool,int);
		void SetTraceCFG();
		
		void SetAlpha(int,double*);

		/* Trace boundary stuff */
		/* Maximum radial distance (default = 1000)*/
		void SetTraceMaxR(double);
		double GetTraceMaxR();


		/* Spheroid/Sphere radii for suface/1 bar level*/
		void SetSurfaceSpheroidR(double,double);
		void GetSurfaceSpheroidR(double*,double*);
		void SetSurfaceSphereR(double);
		double GetSurfaceSphereR();
		void SetSurfaceIsSphere(bool);
		bool GetSurfaceIsSphere();

		/* Spheroid/sphere radii for the ionosphere */
		void SetIonosphereSpheroidR(double,double);
		void GetIonosphereSpheroidR(double*,double*);
		void SetIonosphereSphereR(double);
		double GetIonosphereSphereR();
		void SetIonosphereIsSphere(bool);
		bool GetIonosphereIsSphere();

		/* this is needed to initialize the trace boundaries */
		void SetTraceBoundDefaults();

		/* tracing */
		void TraceField(int*,double**,double**,double**,double**,double**,double**,double**,int**);
		void TraceField();
		void StepVector(double,double,double,double,double*,double*,double*);
		BoolIntTuple ContinueTrace(double,double,double,double*);
		void Step(double,double,double,double*,double*,double*,double*,double*,double*,double*);
		void ReverseElements(int, double*);
		void ReverseElements(int, int*);
		void RKMTrace(	double,double,double,int*,double*,
						double*,double*,double*,double*,double*,double*,int*);
		void FixFootprints(	int,double*,double*,double*,double*,
							double*,double*,double*);
						
		/* get a single field vector */
		void Field(double,double,double,double*,double*,double*);

		/* calculate trace distance,R,Rnorm */
		void CalculateTraceDist(double**);
		void CalculateTraceDist();
		void _CalculateTraceDist();
		void CalculateTraceRnorm(double**);
		void CalculateTraceRnorm();
		void _CalculateTraceRnorm();
	
		/* Calculate footprints */
		void CalculateTraceFP(double**);
		void CalculateTraceFP();
		void _CalculateTraceFP();
		
		/* calculate halpha */
		void CalculateHalpha();
		void CalculateHalpha(double*);
		void CalculateHalpha(double***);
		void CalculateHalpha(double*,double***);

		/* for the conversion between SIII and Mag */
		void SetMagTilt(double);
		double GetMagTilt();
		void SetMagTiltAzimuth(double);
		double GetMagTiltAzimuth();

		/* return things*/
		void GetTraceNstep(int*);
		void GetTrace(double**,double**,double**);
		void GetTrace(double**,double**,double**,double**,double**,double**);
		void GetTraceDist(double**);
		void GetTraceR(double**);
		void GetTraceRnorm(double**);
		void GetTraceFootprints(double**);
		void GetTraceHalpha(double*);	/* python will use this */
		void GetTraceHalpha(double***); /* no idea how to link this to python*/
		
		Trace TracePosition(int,double,double,double);

		/* input coords */
		int n_;
		double *x0_, *y0_, *z0_;  
		int *Date_;
		float *ut_;

		/* trace params */
		int MaxLen_;
		double MaxStep_, MinStep_, InitStep_;
		bool Verbose_;
		int TraceDir_;
		double ErrMax_;
		
		/* this multiplier is basically a hack to make parallel traces
		used to calculate h_alpha trace a bit further into the planet
		so that the interpolation works and there is no extrapolation*/
		double RMultiplier_;

		/* Trace boundaries */
		double MaxR_;
		bool SurfaceIsSphere_, IonosphereIsSphere_;
		double rs_, ri_; //sphere radii
		double as_, bs_, ai_, bi_; //spheroid major/minor axes
		
		/* trace coords */
		int *nstep_;
		double **x_, **y_, **z_;
	
		/* trace fields */
		double **bx_, **by_, **bz_;

		/* trace region - mostly to say whether the trace is:
		* above both surface and ionosphere = 2
		* above surface, below ionosphere = 1
		* below surface, below ionosphere = 0
		* below surface, above ionosphere = -1
		* (that last one might happen if somebody does something odd with the
		* shapes of the surface)
		* */
		int **traceRegion_;

		/* magnetic z-axis tilt (xt) and longitude of tilt (xp)*/
		double xt_;
		double xp_;

		/* trace footprints)*/
		FPstr *fps_;
		FPstr *fpi_;
		EqFPstr *fpe_;
		//double *xfn3_, *yfn3_, *zfn3_;
		//double *xfs3_, *yfs3_, *zfs3_;
		//double *xin3_, *yin3_, *zin3_;
		//double *xis3_, *yis3_, *zis3_;
		//double *xfe3_, *yfe3_, *zfe3_;


		/* trace end points (Magnetic/Dipole coordinates)*/
		//double *xfnm_, *yfnm_, *zfnm_;
		//double *xfsm_, *yfsm_, *zfsm_;
		//double *xinm_, *yinm_, *zinm_;
		//double *xism_, *yism_, *zism_;
		//double *xfem_, *yfem_, *zfem_;


		/* field length, R, Rnorm, Halpha, Footprints */
		int nalpha_;
		double *alpha0_, *alpha1_;
		double Delta_;
		double **S_;
		double **R_;
		double **Rnorm_;
		double *Halpha_;
		double ***Halpha3D_;
		double **FP_;
	
	private:
		/* this is the number of field contributions */
		int nf_;
		std::vector<FieldFuncPtr> Funcs_;

		/* booleans to tell the object what has been done */
		bool inputPos_;
		bool tracedField_,allocTrace_;
		bool hasFootprints_,allocFootprints_,allocFootprintStr_;
		bool allocEqFP_;
		bool hasDist_,allocDist_;
		bool hasRnorm_,allocRnorm_;
		bool hasHalpha_,allocHalpha_, allocHalpha3D_;
		bool allocAlpha_;

		
	
		/* hidden trace functions */
		void _TraceField();

		/* halpha functions */
		bool _CheckHalpha();
		void _CalculateHalpha();
		void _CalculateTraceHalpha(int,int,double*);
		void _CalculateHalphaStartPoints(int i, int j,
							double *xe0, double *ye0, double *ze0,
							double *xe1, double *ye1, double *ze1);
	
};

} /* namespace jupitermag */




#endif
#endif
