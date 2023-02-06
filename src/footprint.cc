#include "footprint.h"

/******************************************************
 *
 *	NAME: footprints(n,x,y,z,a,b,xt,xp,xfn,yfn,zfn,xfs,yfs,zfs)
 *
 *	DESCRIPTION: Attempt to find the footprints on a
 * 	spheroidal surface.
 *
 *
 *
 *
 *
 *
 *
 *******************************************************/
void footprints(int n, double *x, double *y, double *z, 
				double a, double b, double xt, double xp,
				double *xfn, double *yfn, double *zfn,
				double *xfs, double *yfs, double *zfs) {

	/* determine whether each end is a footprint */
	double r,t,rp,tp,rhop,zp;
	
	r = sqrt(x[0]*x[0] + y[0]*y[0] + z[0]*z[0]);
	t = asin(x[0]/r);
	rhop = a*cos(t);
	zp = b*sin(t);
	rp = sqrt(rhop*rhop + zp*zp);
	bool begfp = (r <= rp);
	printf("%f %f\n",r,rp);
	printf("%f %f %f\n",x[0],y[0],z[0]);
	printf("%f %f \n",x[n-1],x[n-2]);
	r = sqrt(x[n-1]*x[n-1] + y[n-1]*y[n-1] + z[n-1]*z[n-1]);
	t = asin(x[n-1]/r);
	rhop = a*cos(t);
	zp = b*sin(t);
	rp = sqrt(rhop*rhop + zp*zp);
	bool endfp = (r <= rp);
	printf("%f %f\n",r,rp);
	/* determine which is north and which is south */
	int begns = 0;
	int endns = 0;
	int indmxr = 0;
	double tx = 0.0;
	double ty = 0.0;
	double tz = 0.0;
	int i;
	printf("Here %d %d\n",begfp, endfp);
	if (begfp && endfp) {
		/* if both ends are footprints - then the north one
		should be the one with the largest positive z-value */
		if (z[0] > z[n-1]) {
			begns = 1;
			endns = -1;
		} else {
			begns = -1;
			endns = 1;
		}

		/* find the maximum r index */
		rp = sqrt(x[0]*x[0] + y[0]*y[0] + z[0]*z[0]);
		for (i=1;i<n;i++) {
			r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
			if (r > rp) {
				rp = r;
				indmxr = i;
			}
		}
	} else {
		if (begfp) {
			/*convert to mag coordinates */
			SIIItoMag(x[0],y[0],z[0],xt,xp,&tx,&ty,&tz);
			if (tz >= 0) {
				begns = 1;
			} else {
				begns = -1;
			}
		}
		if (endfp) {
			/*convert to mag coordinates */
			SIIItoMag(x[n-1],y[n-1],z[n-1],xt,xp,&tx,&ty,&tz);
			if (tz > 0) {
				endns = 1;
			} else {
				endns = -1;
			}
		}
	}


	/* find the first footprint indices */
	int i0, i1;
	double xb, yb, zb;
	if (begfp) {
		i0 = 0;
		if (endfp) {
			i1 = indmxr;
		} else {
			i1 = n - 1;
		}
		findFootprint(x,y,z,i0,i1,a,b,&xb,&yb,&zb);
	} else {
		xb = NAN;
		yb = NAN;
		zb = NAN;
	}

	double xe, ye, ze;
	if (endfp) {
		i0 = n - 1;
		if (begfp) {
			i1 = indmxr;
		} else {
			i1 = 0;
		}
		findFootprint(x,y,z,i0,i1,a,b,&xe,&ye,&ze);
	} else {
		xe = NAN;
		ye = NAN;
		ze = NAN;
	}

	/* now put the beginning/end footprints in the correct 
	north/south fopotprint output*/
	if (begns == 1) {
		*xfn = xb;
		*yfn = yb;
		*zfn = zb;
		*xfs = xe;
		*yfs = ye;
		*zfs = ze;

	} else {
		*xfn = xe;
		*yfn = ye;
		*zfn = ze;
		*xfs = xb;
		*yfs = yb;
		*zfs = zb;
	}

}

double planetRadius(double x, double y, double z, double a, double b) {

	double r, t;
	r = sqrt(x*x + y*y + z*z);
	t = asin(x/r);
	
	return thetaPlanetRadius(t,a,b);

}

double thetaPlanetRadius(double theta, double a, double b) {

	double rhop = a*cos(theta);
	double zp = b*sin(theta);
	return sqrt(rhop*rhop + zp*zp);
}

bool isCrossing(double x0, double y0, double z0,
				double x1, double y1, double z1,
				double a, double b) {
	
	double r0, r1, rp, t0, t1, tmid;

	r0 = sqrt(x0*x0 + y0*y0 + z0*z0);
	r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
	t0 = asin(z0/r0);
	t1 = asin(z1/r1);
	tmid = 0.5*(t0 + t1);
	rp = thetaPlanetRadius(tmid,a,b);

	double rmin = std::min(r0,r1);
	double rmax = std::max(r0,r1);

	return (rp >= rmin) && (rp <= rmax);

}

void interpCrossing(double x0, double y0, double z0,
					double x1, double y1, double z1,
					double a, double b,
					double *xfp, double *yfp, double *zfp) {
	printf("Interp\n");
	double rho0 = sqrt(x0*x0 + y0*y0);
	double rho1 = sqrt(x1*x1 + y1*y1);

	double t0 = atan(z0/rho0);
	double t1 = atan(z1/rho1);

	double deg = 0.01745;
	double tp0, tp1, rhop0, rhop1, zp0, zp1;
	if (fabs(t1 - t0) < deg) {
		tp0 = 0.5*(t0 + t1) + deg/2.0;
		tp1 = 0.5*(t0 + t1) - deg/2.0;
	} else {
		tp0 = t0;
		tp1 = t1;
	}
	rhop0 = a*cos(tp0);
	rhop1 = a*cos(tp1);
	zp0 = b*sin(tp0);
	zp1 = b*sin(tp1);

	double mt = (z1 - z0)/(rho1 - rho0);
	double ct = z0 - mt*rho0;

	double mr = (zp1 - zp0)/(rhop1 - rhop0);
	double cr = zp0 - mr*rhop0;

	double mx = (x1 - x0)/(z1 - z0);
	double cx = x0 - mx*z0;

	double my = (y1 - y0)/(y1 - y0);
	double cy = y0 - my*z0;

	double rhofp = (ct - cr)/(mr - mt);
	*zfp = rhofp*mt + ct;
	*xfp = mx*(*zfp) + cx;
	*yfp = my*(*yfp) + cy;
	

}

void findFootprint(	double *x, double *y, double *z,
					int starti, int endi, 
					double a, double b,
					double *xfp, double *yfp, double *zfp) {
	printf("findFootprint\n");
	/* make sure that there are enough points to scan */
	if (abs(starti-endi) < 1) {
		*xfp = NAN;
		*yfp = NAN;
		*zfp = NAN;
		return;
	}

	/* get the direction */
	int dir;
	if (starti > endi) {
		dir = -1;
	} else {
		dir = 1;		*xfp = NAN;
		*yfp = NAN;
		*zfp = NAN;
		return;
	}

	/* find the indices surrounding the crossing */
	int i, i0 = -1, i1 = -1;
	double r0, r1, rp, t0, t1, tmid;
	r0 = sqrt(x[starti]*x[starti] + y[starti]*y[starti] + z[starti]*z[starti]);
	t0 = asin(z[starti]/r0);
	for (i=starti;i<endi;i+=dir) {
		if (isCrossing(x[i],y[i],z[i],x[i+dir],y[i+dir],z[i+dir],a,b)) {
			i0 = i;
			i1 = i + dir;
			break;
		}
	}

	if (i0 == -1) {
		*xfp = NAN;
		*yfp = NAN;
		*zfp = NAN;
		return;		
	}

	/* interpolate to find the footprint */
	interpCrossing(x[i0],y[i0],z[i0],x[i1],y[i1],z[i1],a,b,xfp,yfp,zfp);

}

void eqfootprints(	int n, double *x, double *y, double *z,
					double *xfe, double *yfe, double *zfe, 
					double *L, double *Lon) {
	
	/* find the furthest point along the field line */
	int indmxr = -1;
	int i;
	double rmx = 0.0;
	double r;
	for (i=0;i<n;i++) {
		r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
		if (r > rmx) {
			indmxr = i;
			rmx = r;
		}
	}

	/* set the footprints */
	*xfe = x[indmxr];
	*yfe = y[indmxr];
	*zfe = z[indmxr];

	/* set L shell and calculate longitude*/
	*L = rmx;
	*Lon = 180*atan2((*yfe),(*xfe))/M_PI;


}