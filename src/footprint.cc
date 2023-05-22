#include "footprint.h"


void FillFPOutputArray(int n, FPstr *fpi, FPstr *fps, EqFPstr *fpe, double **FP) {

	int i;
	for (i=0;i<n;i++) {
		FP[i][0] = fpi[i].xn3;
		FP[i][1] = fpi[i].yn3;
		FP[i][2] = fpi[i].zn3;

		FP[i][3] = fpi[i].xs3;
		FP[i][4] = fpi[i].ys3;
		FP[i][5] = fpi[i].zs3;

		FP[i][6] = fpi[i].xnm;
		FP[i][7] = fpi[i].ynm;
		FP[i][8] = fpi[i].znm;

		FP[i][9] = fpi[i].xsm;
		FP[i][10] = fpi[i].ysm;
		FP[i][11] = fpi[i].zsm;

		FP[i][12] = fps[i].xn3;
		FP[i][13] = fps[i].yn3;
		FP[i][14] = fps[i].zn3;

		FP[i][15] = fps[i].xs3;
		FP[i][16] = fps[i].ys3;
		FP[i][17] = fps[i].zs3;

		FP[i][18] = fps[i].xnm;
		FP[i][19] = fps[i].ynm;
		FP[i][20] = fps[i].znm;

		FP[i][21] = fps[i].xsm;
		FP[i][22] = fps[i].ysm;
		FP[i][23] = fps[i].zsm;

		FP[i][24] = fpe[i].x3;
		FP[i][25] = fpe[i].y3;
		FP[i][26] = fpe[i].z3;

		FP[i][27] = fpe[i].xm;
		FP[i][28] = fpe[i].ym;
		FP[i][29] = fpe[i].zm;

		FP[i][30] = fpi[i].lonn;
		FP[i][31] = fpi[i].latn;

		FP[i][32] = fpi[i].mlonn;
		FP[i][33] = fpi[i].mlatn;

		FP[i][34] = fpi[i].lons;
		FP[i][35] = fpi[i].lats;

		FP[i][36] = fpi[i].mlons;
		FP[i][37] = fpi[i].mlats;

		FP[i][38] = fps[i].lonn;
		FP[i][39] = fps[i].latn;

		FP[i][40] = fps[i].mlonn;
		FP[i][41] = fps[i].mlatn;

		FP[i][42] = fps[i].lons;
		FP[i][43] = fps[i].lats;

		FP[i][44] = fps[i].mlons;
		FP[i][45] = fps[i].mlats;

		FP[i][46] = fpe[i].lshell;
		FP[i][47] = fpe[i].mlone;
		FP[i][48] = fpe[i].fllen;

	}


}

void calculateFootprints(int n, double *x, double *y, double *z,
							double a, double b, double xt, double xp,
							FPstr *fp) {
	
	/* check if each end is a footprint or not*/
	bool begfp = _posisfp(x[0],y[0],z[0],a,b);
	bool endfp = _posisfp(x[n-1],y[n-1],z[n-1],a,b);

	/* work out whether each footprint is north or south */
	int begns, endns, indmxr;
	_nsends(n,x,y,z,xt,xp,begfp,endfp,&begns,&endns,&indmxr);

	/* get the beginning footprint */
	double xb, yb, zb;
	_getbegfp(n,x,y,z,begfp,endfp,indmxr,a,b,&xb,&yb,&zb);

	/* get the ending footprint */
	double xe, ye, ze;
	_getendfp(n,x,y,z,begfp,endfp,indmxr,a,b,&xe,&ye,&ze);

	/* now put the beginning/end footprints in the correct 
	north/south fopotprint output*/
	if (begns == 1) {
		fp->xn3 = xb;
		fp->yn3 = yb;
		fp->zn3 = zb;
		fp->xs3 = xe;
		fp->ys3 = ye;
		fp->zs3 = ze;
	} else {
		fp->xn3 = xe;
		fp->yn3 = ye;
		fp->zn3 = ze;
		fp->xs3 = xb;
		fp->ys3 = yb;
		fp->zs3 = zb;
	}

	/* fill all of the rest of the fields */
	_fillfp(fp,xt,xp);

}

void _fillfp(FPstr *fp, double xt, double xp) {

	/* calcualte mag coords */
	SIIItoMag(fp->xn3,fp->yn3,fp->zn3,xt,xp,&(fp->xnm),&(fp->ynm),&(fp->znm));
	SIIItoMag(fp->xs3,fp->ys3,fp->zs3,xt,xp,&(fp->xsm),&(fp->ysm),&(fp->zsm));

	/* calculate lats and lons */
	double rn, rs;
	rn = sqrt(fp->xn3*fp->xn3 + fp->yn3*fp->yn3 + fp->zn3*fp->zn3);
	rs = sqrt(fp->xs3*fp->xs3 + fp->ys3*fp->ys3 + fp->zs3*fp->zs3);
	
	fp->latn = asin(fp->zn3/rn)*rad2deg;
	fp->lats = asin(fp->zs3/rs)*rad2deg;
	fp->mlatn = asin(fp->znm/rn)*rad2deg;
	fp->mlats = asin(fp->zsm/rs)*rad2deg;

	fp->lonn = atan2(fp->yn3,fp->xn3)*rad2deg;
	fp->lons = atan2(fp->ys3,fp->xs3)*rad2deg;
	fp->mlonn = atan2(fp->ynm,fp->xnm)*rad2deg;
	fp->mlons = atan2(fp->ysm,fp->xsm)*rad2deg;

}

void _getbegfp(	int n, double *x, double *y, double *z,
				bool begfp, bool endfp, int indmxr,
				double a, double b,
				double *xfp, double *yfp, double *zfp) {

	/* get the footprint at the beginning of the array */

	int i0, i1;
	xfp[0] = NAN;
	yfp[0] = NAN;
	zfp[0] = NAN;
	if (begfp) {
		i0 = 0;
		if (endfp) {
			i1 = indmxr;
		} else {
			i1 = n - 1;
		}
		findFootprint(x,y,z,i0,i1,a,b,xfp,yfp,zfp);
	}
}

void _getendfp(	int n, double *x, double *y, double *z,
				bool begfp, bool endfp, int indmxr,
				double a, double b,
				double *xfp, double *yfp, double *zfp) {

	/* get the footprint at the beginning of the array */

	int i0, i1;
	xfp[0] = NAN;
	yfp[0] = NAN;
	zfp[0] = NAN;
	if (endfp) {
		i0 = n - 1;
		if (begfp) {
			i1 = indmxr;
		} else {
			i1 = 0;
		}
		findFootprint(x,y,z,i0,i1,a,b,xfp,yfp,zfp);
	}
}

bool _posisfp(double x, double y, double z, double a, double b) {
	/* determine whether a position (the end of a trace)
	 * is a footprint by whether is is below the surface*/

	double r,t,rp,tp,rhop,zp;
	r = sqrt(x*x + y*y + z*z);
	t = asin(x/r);
	rhop = a*cos(t);
	zp = b*sin(t);
	rp = sqrt(rhop*rhop + zp*zp);
	return r <= rp;
}

int _maxR(double n, double *x, double *y, double *z) {
	/* find the index of the maximum in R*/
	double rmx = -1.0;
	double r;
	int imx = -1;
	int i;
	for (i=0;i<n;i++) {
		r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
		if (r > rmx) {
			rmx = r;
			imx = i;
		}
	}
	return imx;

}

void _nsends(int n, double *x, double *y, double *z,
			double xt, double xp,
			bool begfp, bool endfp,
			int *begns, int *endns, int *indmxr) {
	/* find the maximum r index */
	indmxr[0] = _maxR(n,x,y,z);
	
	/* determine which is north and which is south */
	begns[0] = 0;
	endns[0] = 0;
	
	double tx = 0.0;
	double ty = 0.0;
	double tz = 0.0;
	int i;
	//printf("Here %d %d\n",begfp, endfp);
	if (begfp && endfp) {
		/* if both ends are footprints - then the north one
		should be the one with the largest positive z-value */
		if (z[0] > z[n-1]) {
			begns[0] = 1;
			endns[0] = -1;
		} else {
			begns[0] = -1;
			endns[0] = 1;
		}
	
	} else {
		if (begfp) {
			/*convert to mag coordinates */
			SIIItoMag(x[0],y[0],z[0],xt,xp,&tx,&ty,&tz);
			if (tz >= 0) {
				begns[0] = 1;
			} else {
				begns[0] = -1;
			}
		}
		if (endfp) {
			/*convert to mag coordinates */
			SIIItoMag(x[n-1],y[n-1],z[n-1],xt,xp,&tx,&ty,&tz);
			if (tz > 0) {
				endns[0] = 1;
			} else {
				endns[0] = -1;
			}
		}
	}	 



}

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
	//printf("%f %f\n",r,rp);
	//printf("%f %f %f\n",x[0],y[0],z[0]);
	//printf("%f %f \n",x[n-1],x[n-2]);
	r = sqrt(x[n-1]*x[n-1] + y[n-1]*y[n-1] + z[n-1]*z[n-1]);
	t = asin(x[n-1]/r);
	rhop = a*cos(t);
	zp = b*sin(t);
	rp = sqrt(rhop*rhop + zp*zp);
	bool endfp = (r <= rp);
	//printf("%f %f\n",r,rp);
	/* determine which is north and which is south */
	int begns = 0;
	int endns = 0;
	int indmxr = 0;
	double tx = 0.0;
	double ty = 0.0;
	double tz = 0.0;
	int i;
	//printf("Here %d %d\n",begfp, endfp);
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
	double xb = NAN, yb = NAN, zb = NAN;
	begfp = false;
	if (begfp) {
		i0 = 0;
		if (endfp) {
			i1 = indmxr;
		} else {
			i1 = n - 1;
		}
		findFootprint(x,y,z,i0,i1,a,b,&xb,&yb,&zb);
	}

	double xe = NAN, ye = NAN, ze = NAN;
	endfp = false;
	if (endfp) {
		i0 = n - 1;
		if (begfp) {
			i1 = indmxr;
		} else {
			i1 = 0;
		}
		findFootprint(x,y,z,i0,i1,a,b,&xe,&ye,&ze);
	}
	//printf("%f %f %f %f %f %f\n",xe,ye,ze,xb,yb,zb);
	/* now put the beginning/end footprints in the correct 
	north/south fopotprint output*/
	if (begns == 1) {
		xfn[0] = xb;
		yfn[0] = yb;
		zfn[0] = zb;
		xfs[0] = xe;
		yfs[0] = ye;
		zfs[0] = ze;

	} else {
		xfn[0] = xe;
		yfn[0] = ye;
		zfn[0] = ze;
		xfs[0] = xb;
		yfs[0] = yb;
		zfs[0] = zb;
	}
	//printf("FP: %f %f %f %f %f %f\n",*xfn,*yfn,*zfn,*xfs, *yfs,*zfs);
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
	//printf("%f %f %f %f\n",r0,r1,rp,tmid);
	double rmin = std::min(r0,r1);
	double rmax = std::max(r0,r1);

	return (rp >= rmin) && (rp <= rmax);

}

void interpCrossing(double x0, double y0, double z0,
					double x1, double y1, double z1,
					double a, double b,
					double *xfp, double *yfp, double *zfp) {
	//printf("Interp\n");
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
	*yfp = my*(*zfp) + cy;
	

}

void findFootprint(	double *x, double *y, double *z,
					int starti, int endi, 
					double a, double b,
					double *xfp, double *yfp, double *zfp) {
	//printf("findFootprint\n");
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
		dir = 1;		
	}

	
	/* find the indices surrounding the crossing */
	int i, i0 = -1, i1 = -1;
	double r0, r1, rp, t0, t1, tmid;
	r0 = sqrt(x[starti]*x[starti] + y[starti]*y[starti] + z[starti]*z[starti]);
	t0 = asin(z[starti]/r0);
	if (dir == 1) {
		for (i=starti;i<endi;i++) {
			if (isCrossing(x[i],y[i],z[i],x[i+dir],y[i+dir],z[i+dir],a,b)) {
				i0 = i;
				i1 = i + dir;
				break;
			}
		}
	} else {
		for (i=starti;i>endi;i--) {
			if (isCrossing(x[i],y[i],z[i],x[i-1],y[i-1],z[i-1],a,b)) {
				i0 = i;
				i1 = i + dir;
				break;
			}
		}		
	}
	//printf("%d %d %d %d %d\n",dir,starti,endi,i0,i1);

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

double _fllen(int n, double *x, double *y, double *z) {
	int i;
	double s = 0.0;
	double dx, dy, dz;
	for (i=0;i<n-1;i++) {
		dx = x[i] - x[i+1];
		dy = y[i] - y[i+1];
		dz = z[i] - z[i+1];
		s += sqrt(dx*dx + dy*dy + dz*dz);
	}
	return s;
}

void calculateEquatorialFootprints(int n, double *x, double *y, double *z, 
						double xt, double xp, EqFPstr *efp) {

	/* find the furthest point along the field line*/
	int indmxr = _maxR(n,x,y,z);

	efp->x3 = x[indmxr];
	efp->y3 = y[indmxr];
	efp->z3 = z[indmxr];

	/* convert to mag */
	SIIItoMag(efp->x3,efp->y3,efp->z3,xt,xp,&(efp->xm),&(efp->ym),&(efp->zm));

	/* lshell, lon, field line length */
	efp->fllen = _fllen(n,x,y,z);
	efp->lshell = sqrt(x[indmxr]*x[indmxr] + y[indmxr]*y[indmxr] + z[indmxr]*z[indmxr]);
	efp->mlone = atan2(y[indmxr],x[indmxr]);
}