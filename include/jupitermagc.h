#ifndef __LIBJUPITERMAG_H__
#define __LIBJUPITERMAG_H__
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define LIBJUPITERMAG_VERSION_MAJOR 1
#define LIBJUPITERMAG_VERSION_MINOR 1
#define LIBJUPITERMAG_VERSION_PATCH 0

#define INTERNALFIELD_VERSION_MAJOR 1
#define INTERNALFIELD_VERSION_MINOR 1
#define INTERNALFIELD_VERSION_PATCH 0
#define LIBCON2020_VERSION_MAJOR 0
#define LIBCON2020_VERSION_MINOR 1
#define LIBCON2020_VERSION_PATCH 0
#define LIBSPLINE_VERSION_MAJOR 0
#define LIBSPLINE_VERSION_MINOR 0
#define LIBSPLINE_VERSION_PATCH 1
#define M_PI		3.14159265358979323846

#define __LIBCON2020_H__
#define __LIBINTERNALFIELD_H__
#define __LIBSPLINE_H__
#define deg2rad M_PI/180.0;

	bool TraceField(int n, double *x0, double *y0, double *z0,
					const char *IntFunc, int nExt, char **ExtFunc,
					int MaxLen, double MaxStep, double InitStep,
					double MinStep, double ErrMax, double Delta,
					bool Verbose, int TraceDir,
					int *nstep,
					double **x, double **y, double **z,
					double **Bx, double **By, double **Bz,
					double **R, double **S, double **Rnorm, double **FP,
					int nalpha, double *alpha, double *halpha);
	void ModelField(double p0, double p1, double p2, 
					const char *internal, const char *external, 
					bool CartIn, bool CartOut,
					double *B0, double *B1, double *B2);

	void ModelFieldArray(	int n, double *p0, double *p1, double *p2, 
							const char *internal, const char *external, 
							bool CartIn, bool CartOut,
							double *B0, double *B1, double *B2);
	/* these wrappers can be used to get the magnetic field vectors */
	void Con2020FieldArray(int n, double *p0, double *p1, double *p2,
					double *B0, double *B1, double *B2);
	
	void Con2020Field(double p0, double p1, double p2,
			double *B0, double *B1, double *B2);


	void GetCon2020Params(double *mui, double *irho, double *r0, double *r1,
					double *d, double *xt, double *xp, char *eqtype,
					bool *Edwards, bool *ErrChk, bool *CartIn, bool *CartOut, 
					bool *smooth, double *DeltaRho, double *DeltaZ,
					double *g, char *azfunc, double *wO_open, double *wO_oc,
					double *thetamm, double *dthetamm, double *thetaoc, double *dthetaoc);
						
	
	void SetCon2020Params(double mui, double irho, double r0, double r1,
					double d, double xt, double xp, const char *eqtype,
					bool Edwards, bool ErrChk, bool CartIn, bool CartOut, 
					bool smooth, double DeltaRho, double DeltaZ,
					double g, const char *azfunc, double wO_open, double wO_oc,
					double thetamm, double dthetamm, double thetaoc, double dthetaoc);

	void Con2020AnalyticField(	int n, double a, 
							double *rho, double *z, 
							double *Brho, double *Bz);

	void Con2020AnalyticFieldSmooth(	int n, double a, 
							double *rho, double *z, 
							double *Brho, double *Bz);


/***************************************************************
*
*   NAME : ScalarPotentialSmallRho(rho,z,a,mui2,D)
*
*   DESCRIPTION : Calcualte the small rho approximation
*       of the scalar potential accoring to Edwards et al.,
*       2001 (equation 8).
*
*   INPUTS : 
*       double  rho     Cylindrical rho coordinate (in disc 
*                       coordinate system, Rj)
*       double  z       z-coordinate, Rj
*       double  a       inner edge of semi-infinite current 
*                       sheet, Rj
*       double mui2     mu_0 I_0 /2 parameter (default 139.6 nT)
*       double D        Current sheet half-thickness, Rj
*
***************************************************************/
	double ScalarPotentialSmallRho( double rho, double z, double a,
									double mui2, double D);


/***************************************************************
*
*   NAME : ScalarPotentialLargeRho(rho,z,a,mui2,D,deltaz)
*
*   DESCRIPTION : Calcualte the large rho approximation
*       of the scalar potential accoring to Edwards et al.,
*       2001 (equation 12).
*
*   INPUTS : 
*       double  rho     Cylindrical rho coordinate (in disc 
*                       coordinate system, Rj)
*       double  z       z-coordinate, Rj
*       double  a       inner edge of semi-infinite current 
*                       sheet, Rj
*       double mui2     mu_0 I_0 /2 parameter (default 139.6 nT)
*       double D        Current sheet half-thickness, Rj
*       double deltaz   Scale length over which to smooth 4th
*                        term of the equation
*
***************************************************************/
	double ScalarPotentialLargeRho( double rho, double z, double a,
									double mui2, double D, double deltaz);


/***************************************************************
*
*   NAME : ScalarPotential(rho,z,a,mui2,D,deltarho,deltaz)
*
*   DESCRIPTION : Calculate the small/large rho approximation
*       of the scalar potential accoring to Edwards et al.,
*       2001 (equations 8 & 12).
*
*   INPUTS : 
*       double  rho     Cylindrical rho coordinate (in disc 
*                       coordinate system, Rj)
*       double  z       z-coordinate, Rj
*       double  a       inner edge of semi-infinite current 
*                       sheet, Rj
*       double mui2     mu_0 I_0 /2 parameter (default 139.6 nT)
*       double D        Current sheet half-thickness, Rj
*       double deltarho Sc}
extern "C" {ale length to smoothly transition from
*                       small to large rho approx
*       double deltaz   Scale length over which to smooth 4th
*                        term of the equation
*
***************************************************************/
	double ScalarPotential( double rho, double z, double a,
							double mui2, double D, 
							double deltarho, double deltaz);

	/*************************************************************
	*
	*	NAME: f_theta(thetai)
	*
	*	DESCRIPTION: Equation 5 of Cowley et al., 2008
	*
	*	INPUTS:
	*		double thetai	colatitude of the ionospheric footprint
	*						in radians!
	*
	*	RETURNS:
	*		double f_theta	1 + 0.25*tan^2 thetai
	*
	*************************************************************/
	double f_thetai(double thetai);

	/*************************************************************
	*
	*	NAME: OmegaRatio(thetai,wO_open,wO_om,thetamm,dthetamm,
	*						thetaoc,dthetaoc)
	*
	*	DESCRIPTION: Ratio of the angular velocity mapped to
	*		thetai to the planetary rotation. Equation 15 of 
	*		Cowley et al., 2008.
	*
	*	INPUTS:
	*		double thetai	colatitude of the ionospheric footprint
	*						in radians!
	*		double wO_open	angular velocity ratio of open flux to
	*						planetary spin
	*		double wO_om	angular velocity ratio of outer magnetosphere
	*						to planetary spin
	*		double thetamm	ionospheric footprint latitude of the 
	*						middle magnetosphere (where plasma 
	*						goes from rigid corotation to subcorotation)
	*						in radians.
	*		double dthetamm	width of the middle magnetosphere in radians.
	*		double thetaoc	ionospheric latitude of the open-closed field
	*						line boundary, in radians.
	*		double dthetaoc	width of the open-closed field line boundary,
	*						in radians.
	*
	*	RETURNS:
	*		double wO		Ratio of plasma angular veloctiy to Jupiter
	*						spin.
	*
	*************************************************************/
	double OmegaRatio(	double thetai, double wO_open, double wO_om,
						double thetamm, double dthetamm,
						double thetaoc, double dthetaoc);

	/*************************************************************
	*
	*	NAME: PedersenCurrent(thetai,g,wO_open,wO_om,thetamm,dthetamm,
	*						thetaoc,dthetsoc)
	*
	*	DESCRIPTION: Calculate the Pedersen current which maps to a
	*		given ionospheric latitude using equation 6 of Cowley et
	*		al., 2008.
	*
	*	INPUTS:
	*		double thetai	colatitude of the ionospheric footprint
	*						in radians!
	*		double g		dipole coefficient, nT.
	*		double wO_open	angular velocity ratio of open flux to
	*						planetary spin
	*		double wO_om	angular velocity ratio of outer magnetosphere
	*						to planetary spin
	*		double thetamm	ionospheric footprint latitude of the 
	*						middle magnetosphere (where plasma 
	*						goes from rigid corotation to subcorotation)
	*						in radians.
	*		double dthetamm	width of the middle magnetosphere in radians.
	*		double thetaoc	ionospheric latitude of the open-closed field
	*						line boundary, in radians.
	*		double dthetaoc	width of the open-closed field line boundary,
	*						in radians.
	*	RETURNS:
	*		double Ihp		Ionospheric Pedersen current.
	*
	*************************************************************/
	double PedersenCurrent(	double thetai, double g, 
						double wO_open, double wO_om,
						double thetamm, double dthetamm,
						double thetaoc, double dthetaoc );				

	/*************************************************************
	*
	*	NAME: ThetaIonosphere(r,theta,g,r0,r1,mui2,D,deltarho,deltaz)
	*
	*	DESCRIPTION: Use the flux functions of the CAN model and a 
	*		dipole field to map the current position to a position
	*		on the ionosphere.
	*
	*	INPUTS:
	*		double r		radial coordinate, Rj.
	*		double theta	colatitude, radians.
	*		double g		dipole coefficient, nT.
	*		double r0		Inner edge of the current sheet, Rj.
	*		double r1		Outer edge of the current sheet, Rj.
	*		double mui2		current parameter, nT.
	*		double D		half-thickness of the current sheet, Rj.
	*		double deltarho	scale distance of the smoothing between
	*						inner and outer approximations, Rj.
	*		double deltaz	scale distance to smooth across the
	*						+/-D boundary, Rj.
	*
	*	RETURNS:
	*		double thetai	Ionospheric latitude in radians.
	*
	*
	*************************************************************/
	double ThetaIonosphere(	double r, double theta, double g,
							double r0, double r1,
							double mui2, double D, 
							double deltarho, double deltaz);

	/*************************************************************
	*
	*	NAME: BphiLMIC(r,theta,g,r0,r1,mui2,D,deltarho,deltaz,
	*					wO_open,wO_om,thetamm,dthetamm,
	*					thetaom,dthetaom)
	*
	*	DESCRIPTION: Calculate the azimuthal field using the LMIC 
	*		model.
	*
	*	INPUTS:
	*		double r		radial coordinate, Rj.
	*		double theta	colatitude, radians.
	*		double g		dipole coefficient, nT.
	*		double r0		Inner edge of the current sheet, Rj.
	*		double r1		Outer edge of the current sheet, Rj.
	*		double mui2		current parameter, nT.
	*		double D		half-thickness of the current sheet, Rj.
	*		double deltarho	scale distance of the smoothing between
	*						inner and outer approximations, Rj.
	*		double deltaz	scale distance to smooth across the
	*						+/-D boundary, Rj.
	*		double wO_open	angular velocity ratio of open flux to
	*						planetary spin
	*		double wO_om	angular velocity ratio of outer magnetosphere
	*						to planetary spin
	*		double thetamm	ionospheric footprint latitude of the 
	*						middle magnetosphere (where plasma 
	*						goes from rigid corotation to subcorotation)
	*						in radians.
	*		double dthetamm	width of the middle magnetosphere in radians.
	*		double thetaoc	ionospheric latitude of the open-closed field
	*						line boundary, in radians.
	*		double dthetaoc	width of the open-closed field line boundary,
	*						in radians.
	*
	*	RETURNS:
	*		double Bphi		Azimuthal field, nT.
	*
	*************************************************************/
	double BphiLMIC(double r, double theta, double g,
							double r0, double r1,
							double mui2, double D, 
							double deltarho, double deltaz,
							double wO_open, double wO_om,
							double thetamm, double dthetamm,
							double thetaoc, double dthetaoc );

	/*************************************************************
	*
	*	NAME: BphiIonosphere(thetai,g,wO_open,wO_om,thetamm,dthetamm,
	*					thetaom,dthetaom)
	*
	*	DESCRIPTION: Calculate the ionospheric azimuthal field using the LMIC 
	*		model.
	*
	*	INPUTS:
	*		double thetai	ionospheric colatitude, radians.
	*		double g		dipole coefficient, nT.
	*		double wO_open	angular velocity ratio of open flux to
	*						planetary spin
	*		double wO_om	angular velocity ratio of outer magnetosphere
	*						to planetary spin
	*		double thetamm	ionospheric footprint latitude of the 
	*						middle magnetosphere boundary (where plasma 
	*						goes from rigid corotation to subcorotation)
	*						in radians.
	*		double dthetamm	width of the middle magnetosphere boundary
	*						in radians.
	*		double thetaoc	ionospheric latitude of the open-closed field
	*						line boundary, in radians.
	*		double dthetaoc	width of the open-closed field line boundary,
	*						in radians.
	*
	*	RETURNS:
	*		double Bphi		Azimuthal field, nT.
	*
	*************************************************************/
	double BphiIonosphere( 	double thetai, double g,
							double wO_open, double wO_om,
							double thetamm, double dthetamm,
							double thetaoc, double dthetaoc );
	void spline(int n0, double *x0, double *y0, 
				int n1, double *x1, double *y1);

/* map of strings to direct field model function pointers */
typedef void (*modelFieldPtr)(double,double,double,double*,double*,double*);

/***********************************************************************
 * NAME : getModelFieldPointer(Model)
 *
 * DESCRIPTION : Function to return a pointer to a wrapper function
 * 			which will provide a single field vector at a single 
 * 			position.
 *		
 * INPUTS : 
 *		const char *Model		Model name (use lower case!).
 *
 * RETURNS :
 *		modelFieldPtr *ptr		Pointer to model wrapper.
 *
 **********************************************************************/
	modelFieldPtr getModelFieldPtr(const char *Model);

/* functions to directly call each model for a single Cartesian vector (this will be used for tracing) */

/***********************************************************************
 * NAME : XXXXXField(x,y,z,Bx,By,Bz)
 *
 * DESCRIPTION : Model wrapper functions which can be passed to the 
 * 			tracing code. Replace XXXXXX with the name of the model...
 *		
 * INPUTS : 
 *		double	x			x coordinate in planetary radii.
 *		double	y			y coordinate in planetary radii.
 *		double	z			z coordinate in planetary radii.
 *
 * OUTPUTS :
 *		double	*Bx			x component of the field (nT).
 *		double	*By			y component of the field (nT).
 *		double	*Bz			z component of the field (nT).
 * 
 **********************************************************************/
	void gsfc15evsField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void vip4Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void v117evField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void gsfc15evField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void gsfc13evField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void vipalField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void jpl15evsField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void u17evField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void jrm09Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void o6Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void o4Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void shaField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void p11aField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void jrm33Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void vit4Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void isaacField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void jpl15evField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void spvField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void soiField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void v2Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void cassini3Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void cassini5Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void z3Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void burton2009Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void v1Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void p1184Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void p11asField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void mh2014Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void cain2003Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void langlais2019Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void gao2021Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1935Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf2005Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf2000Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1950Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1960Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1985Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1945Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1965Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1905Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf2010Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf2020Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1910Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1990Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf2015Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1925Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf2025Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1970Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1930Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1920Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1955Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1995Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1900Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1980Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1940Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1975Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void igrf1915Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void nmohField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void gsfco8fullField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void gsfco8Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void thebault2018m3Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void anderson2010qts04Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void uno2009svdField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void anderson2012Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void thebault2018m1Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void anderson2010dts04Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void anderson2010qField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void anderson2010dField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void anderson2010qshaField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void anderson2010dshaField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void ness1975Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void uno2009Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void anderson2010rField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void thebault2018m2Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void ah5Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void gsfcq3fullField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void gsfcq3Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void umohField(double x, double y, double z,
				double *Bx, double *By, double *Bz);

	/* these wrappers can be used to get the magnetic field vectors */

	/***********************************************************************
	 * NAME : InternalField(n,p0,p1,p2,B0,B1,B2)
	 *
	 * DESCRIPTION : Call the model field function. Coordinates depend 
	 * 		on the model  configuration
	 *		
	 * INPUTS : 
	 * 		int		n			Number of array elements
	 *		double	*p0			x or r coordinate in planetary radii.
	 *		double	*p1			y coordinate in planetary radii or theta 
	 * 							in radians.
	 *		double	*p2			z coordinate in planetary radii or phi
	 * 							in radians.
	 *
	 * OUTPUTS :
	 *		double	*B0			x or r component of the field (nT).
	 *		double	*B1			y or theta component of the field (nT).
	 *		double	*B2			z or phi component of the field (nT).
	 * 
	 **********************************************************************/
	void InternalField(int n, double *p0, double *p1, double *p2,
						double *B0, double *B1, double *B2);

	/***********************************************************************
	 * NAME : InternalFieldDeg(n,p0,p1,p2,MaxDeg,B0,B1,B2)
	 *
	 * DESCRIPTION : Call the model field function. Coordinates depend 
	 * 		on the model  configuration
	 *		
	 * INPUTS : 
	 * 		int		n			Number of array elements
	 *		double	*p0			x or r coordinate in planetary radii.
	 *		double	*p1			y coordinate in planetary radii or theta 
	 * 							in radians.
	 *		double	*p2			z coordinate in planetary radii or phi
	 * 							in radians.
	 * 		int 	MaxDeg		Maximum model degree to use.
	 *
	 * OUTPUTS :
	 *		double	*B0			x or r component of the field (nT).
	 *		double	*B1			y or theta component of the field (nT).
	 *		double	*B2			z or phi component of the field (nT).
	 * 
	 **********************************************************************/
	void InternalFieldDeg(int n, double *p0, double *p1, double *p2,
						int MaxDeg, double *B0, double *B1, double *B2);

	/***********************************************************************
	 * NAME : SetInternalCFG(Model,CartIn,CartOut,MaxDeg)
	 *
	 * DESCRIPTION : Configure the current model.
	 *		
	 * INPUTS : 
	 * 		const char *Model		Model name.
	 * 		bool CartIn				Set to True for Cartesian input
	 * 								coordinates or false for polar.
	 * 		bool CartOut			As above, but for the output.
	 * 		int  MaxDeg				Maximum degree used by model
	 * 
	 **********************************************************************/
	void SetInternalCFG(const char *Model, bool CartIn, bool CartOut, int MaxDeg);

	/***********************************************************************
	 * NAME : GetInternalCFG(Model,CartIn,CartOut,MaxDeg)
	 *
	 * DESCRIPTION : Return the current model configuration.
	 *		
	 * OUTPUTS : 
	 * 		char *Model				Model name.
	 * 		bool CartIn				True for Cartesian input
	 * 								coordinates or false for polar.
	 * 		bool CartOut			As above, but for the output.
	 * 		int  MaxDeg				Maximum degree used by model
	 * 
	 **********************************************************************/
	void GetInternalCFG(char *Model, bool *CartIn, bool *CartOut, int *MaxDeg);

	




/* this will be used for all of the model wrapper functions (configure model first) */
typedef void (*FieldFuncPtr)(double,double,double,double*,double*,double*);




/***************************************************************
*
*   NAME : FluxCan(rho,z,r0,r1,mui2,D,deltarho,deltaz)
*
*   DESCRIPTION : Calculate the flux contribution from the 
* 		CAN current sheet (using Edwards et al. 2001 equations).
*
*   INPUTS : 
*       double  rho     Cylindrical rho coordinate (in disc 
*                       coordinate system, Rj)
*       double  z       z-coordinate, Rj
*       double  r0       inner edge of semi-infinite current 
*                       sheet, Rj
*		double 	r1		inner edge of the outer portion of the 
*						current sheet to be subtracted
*       double mui2     mu_0 I_0 /2 parameter (default 139.6 nT)
*       double D        Current sheet half-thickness, Rj
*       double deltarho Scale length to smoothly transition from
*                       small to large rho approx
*       double deltaz   Scale length over which to smooth 4th
*                        term of the equation
*
***************************************************************/
double FluxCan(	double rho,double z, double r0, double r1,
				double mui2, double D, 
				double deltarho, double deltaz) {

	double A0 = ScalarPotential(rho,z,r0,mui2,D,deltarho,deltaz);
	double A1 = ScalarPotential(rho,z,r1,mui2,D,deltarho,deltaz);

	/* according to Edwards et al., 2001 the flux is simply
	rho times the scalar potential */
	double F = rho*(A0 - A1);

	return F;
}

/***************************************************************
*
*   NAME : FluxDip(r,theta,g)
*
*   DESCRIPTION : Calculate the flux cfunction for a dipole
*
*   INPUTS : 
*       double  r 		radial coordinate, Rj
*       double  theta   Colatitude, Rads
*       double  g		Magnetic dipole coefficient, nT
*
***************************************************************/
double FluxDip(double r, double theta, double g) {

	double sint = sin(theta);
	double F = (g*sint*sint)/r;
	return F;
}




#endif
