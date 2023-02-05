#include "coordconv.h"

void MagtoSIII(	double xm, double ym, double zm, 
				double xt, double xp,
				double *x3, double *y3, double *z3) {

	/* angles should be provided in rads */
	double cosxt = cos(xtr);
	double sinxt = sin(xtr);

	double cosxp = cos(xpr - M_PI);
	double sinxp = sin(xpr - M_PI);

	double xtmp = xm*cosxt - zm*sinxt;
	*x3 = xtmp*cosxp - ym*sinxp;
	*y3 = ym*cosxp + xtmp*sinxp;
	*z3 = xm*sinxt + zm*cosxt;
	 

}

void SIIItoMag(	double x3, double y3, double z3, 
				double xt, double xp,
				double *xm, double *ym, double *zm) {

	/* angles should be provided in rads */
	double cosxt = cos(xtr);
	double sinxt = sin(xtr);

	double cosxp = cos(xpr - M_PI);
	double sinxp = sin(xpr - M_PI);

	double xtmp = x3*cosxp + y3*sinxp;
	*xm = xtmp*cosxt + z3*sinxt;
	*ym = y3*cosxp - xtmp*sinxp;
	*zm = x3*sinxt - z3*cosxt;
	 

}