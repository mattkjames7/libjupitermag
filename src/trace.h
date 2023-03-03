#ifndef __TRACE_H__
#define __TRACE_H__
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include "interptraceclosestpos.h"
#include "footprint.h"
#include "coordconv.h"

//const double deg2rad = M_PI/180.0;
//const double rad2deg = 180.0/M_PI;



/* this will be used for all of the model wrapper functions (configure model first) */
typedef void (*FieldFuncPtr)(double,double,double,double*,double*,double*);



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
		void TraceField(int*,double**,double**,double**,double**,double**,double**,double**);
		void TraceField();
		void StepVector(double,double,double,double,double*,double*,double*);
		bool ContinueTrace(double,double,double,double*);
		void Step(double,double,double,double*,double*,double*,double*,double*,double*,double*);
		void ReverseElements(int, double*);
		void RKMTrace(	double,double,double,int*,double*,
						double*,double*,double*,double*,double*,double*);
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
#endif