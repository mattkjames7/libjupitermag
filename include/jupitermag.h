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
	#include <vector>
#else
	#include <math.h>
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include <stdbool.h>
#endif

typedef void (*modelFieldPtr)(double,double,double,double*,double*,double*);


#define LIBJUPITERMAG_VERSION_MAJOR 1
#define LIBJUPITERMAG_VERSION_MINOR 3
#define LIBJUPITERMAG_VERSION_PATCH 0
#ifdef __cplusplus
extern "C" {
#endif
		void spline(int n0, double *x0, double *y0, 
				int n1, double *x1, double *y1);
	
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
*       double deltarho Scale length to smoothly transition from
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
	/* these wrappers can be used to get the magnetic field vectors */
	void Con2020FieldArray(int n, double *p0, double *p1, double *p2,
					double *B0, double *B1, double *B2);
	
	void Con2020Field(double p0, double p1, double p2,
			double *B0, double *B1, double *B2);


	void GetCon2020Params(double *mui, double *irho, double *r0, double *r1,
					double *d, double *xt, double *xp, char *eqtype,
					bool *Edwards, bool *ErrChk, bool *CartIn, bool *CartOut, 
					bool *smooth, double *DeltaRho, double *DeltaZ,
					double *g, char *azfunc, double *wO_open, double *wO_om,
					double *thetamm, double *dthetamm, double *thetaoc, double *dthetaoc);
						
	
	void SetCon2020Params(double mui, double irho, double r0, double r1,
					double d, double xt, double xp, const char *eqtype,
					bool Edwards, bool ErrChk, bool CartIn, bool CartOut, 
					bool smooth, double DeltaRho, double DeltaZ,
					double g, const char *azfunc, double wO_open, double wO_om,
					double thetamm, double dthetamm, double thetaoc, double dthetaoc);

	void Con2020AnalyticField(	int n, double a, 
							double *rho, double *z, 
							double *Brho, double *Bz);

	void Con2020AnalyticFieldSmooth(	int n, double a, 
							double *rho, double *z, 
							double *Brho, double *Bz);

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
	void cassini11Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void p1184Field(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void p11asField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void kivelson2002bField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void kivelson2002aField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void kivelson2002cField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void weber2022dipField(double x, double y, double z,
				double *Bx, double *By, double *Bz);
	void weber2022quadField(double x, double y, double z,
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
					int *nstep,
					double **x, double **y, double **z,
					double **Bx, double **By, double **Bz,
					double **R, double **S, double **Rnorm, double **FP,
					int nalpha, double *alpha, double *halpha);
#ifdef __cplusplus
}

void MagtoSIII(	double xm, double ym, double zm, 
				double xt, double xp,
				double *x3, double *y3, double *z3);

void SIIItoMag(	double x3, double y3, double z3, 
				double xt, double xp,
				double *xm, double *ym, double *zm);


const double deg2rad = M_PI/180.0;
const double rad2deg = 180.0/M_PI;

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









class Spline {
	public:
		Spline(int,double*,double*);
		Spline(const Spline &);
		~Spline();
		void Interpolate(int,double*,double*);
	
		int n_;
		double *a_, *b_, *c_, *d_;
		double *x_, *y_;
		bool del_;
};

	


/* needed this to fix compilation using mingw32 for some reason*/


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
					Spline Sx, Spline Sy, Spline Sz,	/* Splines converting s to a  vector */
					double x, double y, double z,		/* this is the position along the original field line */
					double bx, double by, double bz);
					
					
void OptimizePos(	double x, double y, double z,
					double bx, double by, double bz,
					double s0, 
					Spline Sx, Spline Sy, Spline Sz,
					double *xc, double *yc, double *zc);







double polyeval(double x, double *c, int d);

double pol1eval(double x, double *c, int d);

/***********************************************************************
 * NAME : j0(x)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j0 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		double 	x	position to calculate J0 at.
 * 
 * RETURNS :
 * 		double j	j0 function evaluated at x.
 * 
 * ********************************************************************/
double j0(double x);

/***********************************************************************
 * NAME : j1(x)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j1 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		double 	x	position to calculate J1 at.
 * 
 * RETURNS :
 * 		double j	j1 function evaluated at x.
 * 
 * ********************************************************************/
double j1(double x);

/***********************************************************************
 * NAME : j0(n,x,j)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j0 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		int 	n	Number of elements in x
 * 		double 	*x	position to calculate J0 at.
 * 
 * OUTPUTS :
 * 		double *j	j0 function evaluated at x.
 * 
 * ********************************************************************/
void j0(int n, double *x, double *j);

/***********************************************************************
 * NAME : j1(n,x,j)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j1 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		int 	n	Number of elements in x
 * 		double 	*x	position to calculate J1 at.
 * 
 * OUTPUTS :
 * 		double *j	j1 function evaluated at x.
 * 
 * ********************************************************************/
void j1(int n, double *x, double *j);

/***********************************************************************
 * NAME : j0(n,x,multx,j)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j0 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		int 	n	Number of elements in x
 * 		double 	*x	position to calculate J0(x*multx) at.
 * 		double multx	Constant to multiply x by
 * 
 * OUTPUTS :
 * 		double *j	j0 function evaluated at x*multx.
 * 
 * ********************************************************************/
void j0(int n, double *x, double multx, double *j);

/***********************************************************************
 * NAME : j1(n,x,multx,j)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j1 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		int 	n	Number of elements in x
 * 		double 	*x	position to calculate J1(x*multx) at.
 * 		double multx	Constant to multiply x by
 * 
 * OUTPUTS :
 * 		double *j	j1 function evaluated at x*multx.
 * 
 * ********************************************************************/
void j1(int n, double *x, double multx, double *j);






template <typename T> T clip(T x, T mn, T mx) {
	return std::min(mx,std::max(x,mn));
}





template <typename T> T sgn(T x) {
	return (x > 0) - (x < 0);
}


double trap(int n, double *x, double *y);
double trapc(int n, double dx, double *y);


/***********************************************************************
 * NAME : smoothd(z,dz,d)
 * 
 * DESCRIPTION : Smooth fucntion for crossing the current sheet 
 * (replaces the last bit of equation 12 in Edwards et al 2000).
 * 
 * INPUTS : 
 * 		double z	z-coordinate in dipole coordinate system (Rj)
 * 		double dz	Scale of the transition to use (Rj)
 * 		double d	Half thickness of the current sheet.
 * 
 * RETURNS : 
 * 		double out	Smoothed function across the current sheet.
 * 
 * ********************************************************************/
double smoothd(double z, double dz, double d);


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


/* function pointer for input conversion */
class Con2020; /*this is needed for the pointer below */ 
typedef void (Con2020::*InputConvFunc)(int,double*,double*,double*,
						double*,double*,double*,double*,double*,
						double*,double*,double*,double*);
/* Output conversion */
typedef void (Con2020::*OutputConvFunc)(int,double*,double*,double*,
				double*,double*,double*,double*,
				double*,double*,double*,
				double*,double*,double*);

/* Model function */
typedef void (Con2020::*ModelFunc)(double,double,double,double*,double*,double*);

/* analytical approximation equations */
typedef void (Con2020::*Approx)(double,double,double,double,double,double*,double*);

/* azimuthal function */
typedef void (Con2020::*AzimFunc)(double,double,double,double*);

class Con2020 {
	public:
		/* constructors */
		Con2020();
		Con2020(double,double,double,double,double,double,double,const char*,bool,bool,bool,bool);
	
		/* destructor */
		~Con2020();
		
		/* these functions will be used to set the equations used, if
		 * they need to be changed post-initialisation */
		void SetEdwardsEqs(bool);
		void SetEqType(const char*);
		void SetAzCurrentParameter(double);
		void SetRadCurrentParameter(double);
		void SetR0(double);
		void SetR1(double);
		void SetCSHalfThickness(double);
		void SetCSTilt(double);
		void SetCSTiltAzimuth(double);
		void SetErrCheck(bool);
		void SetCartIn(bool);
		void SetCartOut(bool);
		void SetSmooth(bool);
		void SetDeltaRho(double);
		void SetDeltaZ(double);
		void SetOmegaOpen(double);
		void SetOmegaOM(double);
		void SetThetaMM(double);
		void SetdThetaMM(double);
		void SetThetaOC(double);
		void SetdThetaOC(double);
		void SetG(double);
		void SetAzimuthalFunc(const char*);
		
		/* these mamber functions will be the "getter" version of the
		 * above setters */
		bool GetEdwardsEqs();
		void GetEqType(char*);
		double GetAzCurrentParameter();
		double GetRadCurrentParameter();
		double GetR0();
		double GetR1();
		double GetCSHalfThickness();
		double GetCSTilt();
		double GetCSTiltAzimuth();
		bool GetErrCheck();
		bool GetCartIn();
		bool GetCartOut();
		bool GetSmooth();
		double GetDeltaRho();
		double GetDeltaZ();
		double GetOmegaOpen();
		double GetOmegaOM();
		double GetThetaMM();
		double GetdThetaMM();
		double GetThetaOC();
		double GetdThetaOC();
		double GetG();
		void GetAzimuthalFunc(char *);

		/* This function will be used to call the model, it is overloaded
		 * so that we have one for arrays, one for scalars */
		void Field(int,double*,double*,double*,double*,double*,double*);
		void Field(double,double,double,double*,double*,double*);

		/* a function for testing purposes...*/
		void AnalyticField(double,double,double,double*,double*);
		void AnalyticFieldSmooth(double,double,double,double*,double*);

		/* expose the model function pointer*/
		ModelFunc _Model;

	private:
		/* model parameters */
		double mui_,irho_,r0_,r1_,d_,xt_,xp_,disctilt_,discshift_;
		double r0sq_, r1sq_;
		double cosxp_,sinxp_,cosxt_,sinxt_;
		char eqtype_[9];
		bool Edwards_, ErrChk_;
		bool CartIn_,CartOut_;
		double deltaz_,deltarho_;
		bool smooth_;

		/* LMIC parameters*/
		double wO_open_, wO_om_, thetamm_, dthetamm_, thetaoc_, dthetaoc_, g_;
		char azfunc_[10];
		
		/* Bessel function arrays - arrays prefixed with r and z are
		 * to be used for integrals which calcualte Brho and Bz,
		 * respectively */
		int *rnbes_;			/* number of elements for each Bessel function (rho)*/
		int *znbes_;			/* same as above for z */
		double **rlambda_;/* Lambda array to integrate over rho*/
		double **zlambda_;/* Lambda array to integrate over z*/
		double **rj0_lambda_r0_; /* j0(lambda*r0) */
		double **rj1_lambda_rho_;/* j1(lambda*rho) */
		double **zj0_lambda_r0_; /* j0(lambda*r0) */
		double **zj0_lambda_rho_;/* j0(lambda*rho) */

		
		/* arrays to multiply be stuff to be integrated */
		/* these arrays will store the parts of equations 14, 15, 17 
		 * and 18 of Connerny 1981 which only need to be calculated once*/
		double **Eq14_;		/* j0(lambda*r0)*sinh(lamba*d)/lambda */
		double **Eq15_;     /* j0(lambda*r0)*sinh(lamba*d)/lambda */
		double **Eq17_;     /* j0(lambda*r0)*exp(-lamba*d)/lambda */
		double **Eq18_;     /* j0(lambda*r0)/lambda */
		double **ExpLambdaD_;
		

		/* integration step sizes */
		static constexpr double dlambda_ = 1e-4;
		static constexpr double dlambda_brho_ = 1e-4;
		static constexpr double dlambda_bz_ = 5e-5;
		
		/* Arrays containing maximum lambda values */
		double rlmx_array_[6];
		double zlmx_array_[6];


		/* coordinate conversions for positions */
		InputConvFunc _ConvInput;
		void _SysIII2Mag(int,double*,double*,double*,
						double*,double*,double*,double*,double*,
						double*,double*,double*,double*);
		void _PolSysIII2Mag(int,double*,double*,double*,
						double*,double*,double*,double*,double*,
						double*,double*,double*,double*);
		
		
		/* coordinate conversion for magnetic field vector */
		OutputConvFunc _ConvOutput;
		void _BMag2SysIII(int,double*,double*,double*,
							double*,double*,double*,double*,
							double*,double*,double*,
							double*,double*,double*);
		void _BMag2PolSysIII(int,double*,double*,double*,
							double*,double*,double*,double*,
							double*,double*,double*,
							double*,double*,double*);	

		/* Functions to update function pointers */
		void _SetIOFunctions();
		void _SetModelFunctions();
		
							
		/* Azimuthal field */
		AzimFunc _AzimuthalField;
		void _BphiConnerney(int,double*,double*,double*,double*);
		void _BphiConnerney(double,double,double,double*);
		void _BphiLMIC(double,double,double,double*);
		Approx _LargeRho;
		Approx _SmallRho;		
		/* analytic equations */
		void _Analytic(double,double,double,double*,double*,double*);
		void _AnalyticSmooth(double,double,double,double*,double*,double*);
		void _SolveAnalytic(int,double*,double*,double,double*,double*);
		void _AnalyticInner(double,double,double*,double*);
		void _AnalyticOuter(double,double,double*,double*);
		void _AnalyticInnerSmooth(double,double,double*,double*);
		void _AnalyticOuterSmooth(double,double,double*,double*);
		void _LargeRhoConnerney(double,double,double,double,double,double*,double*);
		void _SmallRhoConnerney(double,double,double,double,double,double*,double*);
		void _LargeRhoEdwards(double,double,double,double,double,double*,double*);
		void _LargeRhoEdwardsSmooth(double,double,double,double,double,double*,double*);
		void _SmallRhoEdwards(double,double,double,double,double,double*,double*);
		
		/* integral-related functions */
		void _Integral(double,double,double,double*,double*,double*);
		void _IntegralInner(double, double, double,	double*, double*);
		void _InitIntegrals();
		void _RecalcIntegrals();
		void _DeleteIntegrals();
		void _IntegralChecks(int,double*,int*,int[]);
		void _IntegralCheck(double,int*);
		void _SolveIntegral(int,double*,double*,double*,double*,double*);
		void _IntegrateEq14(int,double,double,double,double*);
		void _IntegrateEq15(int,double,double,double*);
		void _IntegrateEq17(int,double,double,double*);
		void _IntegrateEq18(int,double,double,double*);

		/* hybrid */
		void _Hybrid(double,double,double,double*,double*,double*);
};



/* we want to initialize the model objects with its parameters */
extern Con2020 con2020;


	



/* this is used in both C and C++*/



/* structure for storing the coefficients in memory (replaces binary stuff) */
typedef struct coeffStruct {
    const int len;    
    const int nmax;
    const int ndef;
    const double rscale;
    const int *n;
    const int *m;
    const double *g;
    const double *h;
} coeffStruct;

typedef coeffStruct& (*coeffStructFunc)();


/* list of model names */
std::vector<std::string> getModelNames();

/* model coefficient arrays */
extern coeffStruct& _model_coeff_gsfc15evs();
extern coeffStruct& _model_coeff_vip4();
extern coeffStruct& _model_coeff_v117ev();
extern coeffStruct& _model_coeff_gsfc15ev();
extern coeffStruct& _model_coeff_gsfc13ev();
extern coeffStruct& _model_coeff_vipal();
extern coeffStruct& _model_coeff_jpl15evs();
extern coeffStruct& _model_coeff_u17ev();
extern coeffStruct& _model_coeff_jrm09();
extern coeffStruct& _model_coeff_o6();
extern coeffStruct& _model_coeff_o4();
extern coeffStruct& _model_coeff_sha();
extern coeffStruct& _model_coeff_p11a();
extern coeffStruct& _model_coeff_jrm33();
extern coeffStruct& _model_coeff_vit4();
extern coeffStruct& _model_coeff_isaac();
extern coeffStruct& _model_coeff_jpl15ev();
extern coeffStruct& _model_coeff_spv();
extern coeffStruct& _model_coeff_soi();
extern coeffStruct& _model_coeff_v2();
extern coeffStruct& _model_coeff_cassini3();
extern coeffStruct& _model_coeff_cassini5();
extern coeffStruct& _model_coeff_z3();
extern coeffStruct& _model_coeff_burton2009();
extern coeffStruct& _model_coeff_v1();
extern coeffStruct& _model_coeff_cassini11();
extern coeffStruct& _model_coeff_p1184();
extern coeffStruct& _model_coeff_p11as();
extern coeffStruct& _model_coeff_kivelson2002b();
extern coeffStruct& _model_coeff_kivelson2002a();
extern coeffStruct& _model_coeff_kivelson2002c();
extern coeffStruct& _model_coeff_weber2022dip();
extern coeffStruct& _model_coeff_weber2022quad();
extern coeffStruct& _model_coeff_mh2014();
extern coeffStruct& _model_coeff_cain2003();
extern coeffStruct& _model_coeff_langlais2019();
extern coeffStruct& _model_coeff_gao2021();
extern coeffStruct& _model_coeff_igrf1935();
extern coeffStruct& _model_coeff_igrf2005();
extern coeffStruct& _model_coeff_igrf2000();
extern coeffStruct& _model_coeff_igrf1950();
extern coeffStruct& _model_coeff_igrf1960();
extern coeffStruct& _model_coeff_igrf1985();
extern coeffStruct& _model_coeff_igrf1945();
extern coeffStruct& _model_coeff_igrf1965();
extern coeffStruct& _model_coeff_igrf1905();
extern coeffStruct& _model_coeff_igrf2010();
extern coeffStruct& _model_coeff_igrf2020();
extern coeffStruct& _model_coeff_igrf1910();
extern coeffStruct& _model_coeff_igrf1990();
extern coeffStruct& _model_coeff_igrf2015();
extern coeffStruct& _model_coeff_igrf1925();
extern coeffStruct& _model_coeff_igrf2025();
extern coeffStruct& _model_coeff_igrf1970();
extern coeffStruct& _model_coeff_igrf1930();
extern coeffStruct& _model_coeff_igrf1920();
extern coeffStruct& _model_coeff_igrf1955();
extern coeffStruct& _model_coeff_igrf1995();
extern coeffStruct& _model_coeff_igrf1900();
extern coeffStruct& _model_coeff_igrf1980();
extern coeffStruct& _model_coeff_igrf1940();
extern coeffStruct& _model_coeff_igrf1975();
extern coeffStruct& _model_coeff_igrf1915();
extern coeffStruct& _model_coeff_nmoh();
extern coeffStruct& _model_coeff_gsfco8full();
extern coeffStruct& _model_coeff_gsfco8();
extern coeffStruct& _model_coeff_thebault2018m3();
extern coeffStruct& _model_coeff_anderson2010qts04();
extern coeffStruct& _model_coeff_uno2009svd();
extern coeffStruct& _model_coeff_anderson2012();
extern coeffStruct& _model_coeff_thebault2018m1();
extern coeffStruct& _model_coeff_anderson2010dts04();
extern coeffStruct& _model_coeff_anderson2010q();
extern coeffStruct& _model_coeff_anderson2010d();
extern coeffStruct& _model_coeff_anderson2010qsha();
extern coeffStruct& _model_coeff_anderson2010dsha();
extern coeffStruct& _model_coeff_ness1975();
extern coeffStruct& _model_coeff_uno2009();
extern coeffStruct& _model_coeff_anderson2010r();
extern coeffStruct& _model_coeff_thebault2018m2();
extern coeffStruct& _model_coeff_ah5();
extern coeffStruct& _model_coeff_gsfcq3full();
extern coeffStruct& _model_coeff_gsfcq3();
extern coeffStruct& _model_coeff_umoh();

/* map model names to the structure containing the coefficients */
std::map<std::string,coeffStructFunc> getCoeffMap();

/***********************************************************************
 * NAME : getModelCoeffStruct(Model)
 *
 * DESCRIPTION : Function to return a structure containing model 
        coefficients.
 *		
 * INPUTS : 
 *		std::string Model	Model name (use lower case!).
 *
 * RETURNS :
 *		coeffStructFunc	cstr    Model coefficient function.
 *
 **********************************************************************/
coeffStructFunc getModelCoeffStruct(std::string Model);

/***********************************************************************
 * NAME : getModelCoeffStruct(Model)
 *
 * DESCRIPTION : Function to return a structure containing model 
        coefficients.
 *		
 * INPUTS : 
 *		const char *Model	Model name (use lower case!).
 *
 * RETURNS :
 *		coeffStructFunc	cstr    Model coefficient function.
 *
 **********************************************************************/
coeffStructFunc getModelCoeffStruct(const char *Model);




/* This structure will store the Schmidt coefficients */
struct schmidtcoeffs {
	int n;
	int m;
	double g;
	double h;
};


/***********************************************************************
 * NAME : class Internal
 * 
 * DESCRIPTION : 
 * 		Class which will store the g and h spherical harmonic 
 * 		coefficients for a given model. To obtain the magnetic field,
 * 		use the Field() and FieldCart() member functions.
 * 
 * ********************************************************************/
class Internal {
	public:
		Internal(unsigned char *);
		Internal(const char *);
		Internal(const Internal&);
		~Internal();
	
		/*these four functions will calculate the field in spherical
		 * polar RH system III coordinates.*/
		void Field(int,double*,double*,double*,double*,double*,double*);
		void Field(int,double*,double*,double*,int,double*,double*,double*);
		void Field(double,double,double,double*,double*,double*);
		void Field(double,double,double,int,double*,double*,double*);
		
		/* these will be Cartesian */
		void FieldCart(double,double,double,double*,double*,double*);
		void FieldCart(double,double,double,int,double*,double*,double*);

		/* set current degree */
		void SetDegree(int n);
		int GetDegree();

		
	private:
		/*Schmidt coefficients */
		struct schmidtcoeffs *schc_;
		int nschc_;
		double **Snm_;
		
		/* maximum, default and current degree */
		int nmax_;
		int ndef_;
		int *ncur_;
		
		/* these ones will have Snm_ already multiplied */
		double **g_;
		double **h_;
		
		/* Legendre Polynomial and derivative arrays */
		double **Pnm_, **dPnm_;
		
		/* cosmp and sinmp arrays */
		double *cosmp_, *sinmp_;		
		
		
		/* hack to scale r or x,y,z because some models use a different
		 * definition for the planetary radius - notably the different 
		 * Jupiter models - this should be rpgood/rpbad, where rpgood
		 * is the accepted planetary radius and rpbad is the erroneous
		 * one - this will be then multiplied by r: rnew = r*rscale_
		 * where rscale_ = rgood/rbad */
		double rscale_;
		
		/* functions for initializing the object */
		void _LoadSchmidt(unsigned char*);
		void _LoadSchmidt(coeffStruct );
		void _Schmidt();
		void _CoeffGrids();

		/* This function will calculate the Legendre polynomials */
		void _Legendre(int,double*,double*,int,double***,double***);
		void _Legendre(double,double,int,double**,double**);
		
		/* this function will calculate the magnetic field components in
		 * spherical polar coordinates */
		void _SphHarm(int,double*,double*,double*,double*,double*,double*);
		/* could do with writing a scalar version of this for extra speed */
		void _SphHarm(double,double,double,double*,double*,double*);
		
		void _Cart2Pol(double,double,double,double*,double*,double*);
		void _BPol2BCart(double,double,double,double,double,double*,double*,double*);
		
		bool copy;

		/* initialization */
		bool useptr_;
		bool *init_;
		unsigned char *modelptr_;
		coeffStruct *modelstr_;
		void _Init();
		void _CheckInit();
	
};



/* models! */
extern Internal& gsfc15evs();
extern Internal& vip4();
extern Internal& v117ev();
extern Internal& gsfc15ev();
extern Internal& gsfc13ev();
extern Internal& vipal();
extern Internal& jpl15evs();
extern Internal& u17ev();
extern Internal& jrm09();
extern Internal& o6();
extern Internal& o4();
extern Internal& sha();
extern Internal& p11a();
extern Internal& jrm33();
extern Internal& vit4();
extern Internal& isaac();
extern Internal& jpl15ev();
extern Internal& spv();
extern Internal& soi();
extern Internal& v2();
extern Internal& cassini3();
extern Internal& cassini5();
extern Internal& z3();
extern Internal& burton2009();
extern Internal& v1();
extern Internal& cassini11();
extern Internal& p1184();
extern Internal& p11as();
extern Internal& kivelson2002b();
extern Internal& kivelson2002a();
extern Internal& kivelson2002c();
extern Internal& weber2022dip();
extern Internal& weber2022quad();
extern Internal& mh2014();
extern Internal& cain2003();
extern Internal& langlais2019();
extern Internal& gao2021();
extern Internal& igrf1935();
extern Internal& igrf2005();
extern Internal& igrf2000();
extern Internal& igrf1950();
extern Internal& igrf1960();
extern Internal& igrf1985();
extern Internal& igrf1945();
extern Internal& igrf1965();
extern Internal& igrf1905();
extern Internal& igrf2010();
extern Internal& igrf2020();
extern Internal& igrf1910();
extern Internal& igrf1990();
extern Internal& igrf2015();
extern Internal& igrf1925();
extern Internal& igrf2025();
extern Internal& igrf1970();
extern Internal& igrf1930();
extern Internal& igrf1920();
extern Internal& igrf1955();
extern Internal& igrf1995();
extern Internal& igrf1900();
extern Internal& igrf1980();
extern Internal& igrf1940();
extern Internal& igrf1975();
extern Internal& igrf1915();
extern Internal& nmoh();
extern Internal& gsfco8full();
extern Internal& gsfco8();
extern Internal& thebault2018m3();
extern Internal& anderson2010qts04();
extern Internal& uno2009svd();
extern Internal& anderson2012();
extern Internal& thebault2018m1();
extern Internal& anderson2010dts04();
extern Internal& anderson2010q();
extern Internal& anderson2010d();
extern Internal& anderson2010qsha();
extern Internal& anderson2010dsha();
extern Internal& ness1975();
extern Internal& uno2009();
extern Internal& anderson2010r();
extern Internal& thebault2018m2();
extern Internal& ah5();
extern Internal& gsfcq3full();
extern Internal& gsfcq3();
extern Internal& umoh();


/* map the model names to their model object pointers */
typedef Internal& (*InternalFunc)();
std::map<std::string,InternalFunc> getModelPtrMap();

/* functions to return the pointer to a model object given a string */

/***********************************************************************
 * NAME : getModelObjPointer(Model)
 *
 * DESCRIPTION : Function to return a pointer to a model object.
 *		
 * INPUTS : 
 *		std::string Model	Model name (use lower case!).
 *
 * RETURNS :
 *		InternalFunc ptr		Function pointer to model object.
 *
 **********************************************************************/
InternalFunc getModelObjPointer(std::string Model);

/***********************************************************************
 * NAME : getModelObjPointer(Model)
 *
 * DESCRIPTION : Function to return a pointer to a model object.
 *		
 * INPUTS : 
 *		const char *Model	Model name (use lower case!).
 *
 * RETURNS :
 *		InternalFunc ptr		Function pointer to model object.
 *
 **********************************************************************/
InternalFunc getModelObjPointer(const char *Model);

/* a function to return a list of the models available */
/***********************************************************************
 * NAME : listAvailableModels()
 *
 * DESCRIPTION : Function to return a list of model names available.
 *		
 * RETURNS :
 *		vector<string> Models	Model list.
 *
 **********************************************************************/
std::vector<std::string> listAvailableModels();

/* map of strings to direct field model function pointers */
std::map<std::string,modelFieldPtr> getModelFieldPtrMap();

/* functions to return pointer to model field function */

/***********************************************************************
 * NAME : getModelFieldPointer(Model)
 *
 * DESCRIPTION : Function to return a pointer to a wrapper function
 * 			which will provide a single field vector at a single 
 * 			position.
 *		
 * INPUTS : 
 *		std::string Model		Model name (use lower case!).
 *
 * RETURNS :
 *		modelFieldPtr *ptr		Pointer to model wrapper.
 *
 **********************************************************************/
modelFieldPtr getModelFieldPtr(std::string Model);



/* based upon https://www.lonecpluspluscoder.com/2015/08/13/an-elegant-way-to-extract-keys-from-a-c-map/ */
/***********************************************************************
 * NAME : vector<> listMapKeys(inmap)
 * 
 * DESCRIPTION : List the keys used for a std::map.
 * 
 * INPUTS : 
 * 		map		inmap		std::map instance			
 * 
 * 
 * RETURNS :
 * 		vector	keys		vector object containing a list of the map 
 * 							keys
 * 
 * 
 * 
 * ********************************************************************/
template <typename Tkey, typename Tval> 
std::vector<Tkey> listMapKeys(std::map<Tkey,Tval> const &inmap) {
	std::vector<Tkey> keys;
	for (auto const& element: inmap) {
		keys.push_back(element.first);
	}
	return keys;
}	



/***********************************************************************
 * NAME : class InternalModel
 * 
 * DESCRIPTION : 
 * 		Class which can access all instances of Internal objects.
 * 
 * ********************************************************************/
class InternalModel {
	
	public:
		/* constructor */
		InternalModel();
		
		/* copy constructor */
		InternalModel(const InternalModel&);
		
		/* destructor */
		~InternalModel();
		
		/* Init this function - I would like to remove this if at all possible*/
		void CheckInit();
		void Init();
		
		/* set model parameters */
		void SetCartIn(bool);
		void SetCartOut(bool);
		bool GetCartIn();
		bool GetCartOut();
		void SetModel(const char *);
		void GetModel(char *);
		void SetDegree(int n);
		int GetDegree();

		/* Field functions */
		void Field(int,double*,double*,double*,int,double*,double*,double*);
		void Field(int,double*,double*,double*,double*,double*,double*);
		void Field(double,double,double,int,double*,double*,double*);
		void Field(double,double,double,double*,double*,double*);
				
		/* these objects are the models to use */
		std::map<std::string,Internal*> Models_;
		std::vector<std::string> ModelNames_;

	private:
		Internal *CurrentModel_;
		std::string *CurrentModelName_;


		/* coordinate/field vector rotation */
		bool copy_;
		bool *init_;
		bool *CartIn_;
		bool *CartOut_;
		void _Cart2Pol(int,double*,double*,double*,double*,double*,double*);
		void _Cart2Pol(double,double,double,double*,double*,double*);
		void _BPol2BCart(int,double*,double*,double*,double*,double*,double*,double*,double*);
		void _BPol2BCart(double,double,double,double,double,double*,double*,double*);
};



/* we want to initialize the model objects witht heir parameters */
InternalModel getInternalModel();



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
#endif
	