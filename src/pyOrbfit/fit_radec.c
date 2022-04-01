/* 	$Id: fit_radec.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: fit_radec.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* Fit orbit to RA/DEC file.  First
 * preliminary fit, THEN do a full Marquandt-Levenberg fit.
 * Use sigma matrix of full fit to predict positions and
 * errors using the full nonlinear formulae.  
 * usage:  orbtest3 <datafile> <# params>
 *  datafile contains list of observations (e.g. output of orbtest1)
 *  # params is either 5 (don't fit to gdot) or 6 (fit all params)
 * 8/9/99 gmb
 */
#include "orbfit.h"

char   *help[] = {
    " fit_radec:  Fit KBO orbit to observed astrometric RA/DEC",
    " usage:  fit_radec [-m mpc_error] [-j JPL_file] [-o observatory_file] [-v]",
    "  mpc_error  is arcsec uncertainty to apply to MPC-format",
    "             observations.  Default is 0.2",
    "  JPL_file   is binary ephemeris file.  Default is binEphem.423, or",
    "             a file specified by environment variable ORBIT_EPHEMERIS",
    "  observatory_file   is file with observatory site data.  Default is",
    "             observatories.dat, or a file specified by environment",
    "             variable ORBIT_OBSERVATORIES",
    "  stdin      contains one observation per line, with format",
    "             JD  RA Dec error obscode",
    "  stdout     is a file containing best-fit alpha, beta, etc.",
    "             and its covariance matrix.",
    "  Residuals are dumped to stderr.",
    0
};

void
print_help(char *h[])
{
  int   i;
  for (i = 0; help[i] != 0; i++)
    fprintf (stderr, "%s\n", help[i]);
  exit (1);
}

int
fit_observations(OBSERVATION obsarray[],
		 int nobs,
		 PBASIS *p,
		 double **covar,
		 double *chisq,
		 int *dof,
		 FILE *logfile);

OBSERVATION obsarray[MAXOBS];
int	nobs;

/* Ugly hack to return the covariance matricx in a form that Python can comprehend */
void flatten_cov(double **cov, int ndim, double *cov1d) {
    int i;
    int j;
    for (i=0; i<ndim; i++) {
	for (j=0; j<ndim; j++) {
	    cov1d[ndim*i + j] = cov[i+1][j+1];
	}
    }   
    
}

/* YAUH (yet another ugly hack), this goes in the other direction */
void unflatten_cov(double *cov1d, int ndim, double **cov) {
    int i;
    int j;
    for (i=0; i<ndim; i++) {
	for (j=0; j<ndim; j++) {
	  cov[i+1][j+1] = cov1d[ndim*i + j];
	  printf("i, j, cov[i+1][j+1]: %i, %i, %f\n", i , j,  cov[i+1][j+1]);
	}
    }
    
}

/* Helper function used by Python wrapper */
void add_to_obsarray(OBSERVATION obsarray[], int iobs, OBSERVATION obs) {

    double elat, elon;
    eq_to_ec(obs.thetax,obs.thetay,&elat,&elon,NULL);

    if (iobs==0) {
      double xec, yec, zec;
      /* Use first observation to set the reference frame */
      jd0 = obs.obstime;
      lat0 = elat;
      lon0 = elon;
      
      /* Find location of SSBARY wrt observatory at zero time */
      earth_ssbary(jd0, obs.obscode, &xBary, &yBary, &zBary);

      /* Negate the vector to make it earth->SSBARY*/
      /* And rotate equatorial into the tangent-point coords */
      xBary *= -1.;  yBary *= -1.;  zBary *= -1.;
      xyz_eq_to_ec(xBary, yBary, zBary, &xec, &yec, &zec,NULL);
      xyz_ec_to_proj(xec, yec, zec, &xBary, &yBary, &zBary, lat0, lon0, NULL);
    }
    /* Set time to years after jd0, rotate to tangent plane coords */
    obs.obstime = (obs.obstime-jd0)*DAY;
    ec_to_proj(elat,elon,&obs.thetax,&obs.thetay,lat0,lon0,NULL);
    /* Calculate the position of Earth at this time to avoid doing
     * it many times later: */
    earth3d(obs.obstime, obs.obscode,
	    &obs.xe,&obs.ye, &obs.ze);
    obsarray[iobs]=obs;		
}

void
fit_radec(char *fname, int *nobservations, double *chisqfit, int *ndof, PBASIS *pbasis, ORBIT *orb)
{
  PBASIS p;
  ORBIT  orbit;
  XVBASIS xv;

  double **covar;
  double chisq;
  int i, dof;

  covar = dmatrix(1,6,1,6);
  
  if (read_radec(obsarray, fname, &nobs)) {
/*    fprintf(stderr, "Error reading input observations\n"); */
    exit(1);
  }
  
#include <time.h>
    time_t timettt;
    time(&timettt);
    /* note that ctime returns string with newline at end */
/*    printf("\n#---%s",ctime(&timettt)); */

  printf("# Fitting %d observations\n",nobs); 

  /* Call subroutine to do the actual fitting: */
  fit_observations(obsarray, nobs, &p, covar, &chisq, &dof,stdout);

  *nobservations = nobs;
  *chisqfit = chisq;
  *ndof = dof;
  *pbasis = p;
    

  printf("# Chi-squared of fit: %.2f DOF: %d\n",chisq,dof);
  printf("# Exact a, adot, b, bdot, g, gdot:\n");
  printf("%11.8f %11.8f %11.8f %11.8f %11.8f %11.8f\n",p.a,p.adot,p.b,
	p.bdot, p.g, p.gdot);
  
  pbasis_to_bary(&p, &xv, NULL);

  orbitElements(&xv, &orbit);
  *orb = orbit; 
  printf("# a=%lf AU,e=%lf,i=%lf deg\n",orbit.a, orbit.e, orbit.i);
  {
    double d, dd;
    d = sqrt(xBary*xBary + yBary*yBary + pow(zBary-1/p.g,2.));
    dd = d*d*sqrt(covar[5][5]);
    printf("# Barycentric distance %.3f+-%.3f\n",d,dd); 
  }

  /* Print the covariance matrix to stdout */

  printf("# abg covariance matrix: \n");
  print_matrix(stdout,covar,6,6);


  /* Print out information on the coordinate system */
  printf("#     lat0       lon0       xBary     yBary      zBary   JD0\n");
  printf("%12.7f %12.7f %10.7f %10.7f %10.7f  %.6f\n",
	 lat0/DTOR,lon0/DTOR,xBary,yBary,zBary,jd0);

  /* Dump residuals to stderr */
/*  fprintf(stderr,"Best fit orbit gives:\n");
  fprintf(stderr,"obs  time        x      x_resid       y   y_resid\n"); */
  for (i=0; i<nobs; i++) {
    double x,y;
    kbo2d(&p, &obsarray[i], &x, NULL, &y, NULL);
/*
    fprintf(stderr,"%3d %9.4f %10.3f %7.3f %10.3f %7.3f\n",
	    i, obsarray[i].obstime,
	    obsarray[i].thetax/ARCSEC, (obsarray[i].thetax-x)/ARCSEC,
	    obsarray[i].thetay/ARCSEC, (obsarray[i].thetay-y)/ARCSEC);
*/
  }
  
/* Convert from abg to aei basis */
/* (Code more or less cribbed from abg_to_aei.c) */
  double **covar_xyz = dmatrix(1,6,1,6);
  double **covar_aei = dmatrix(1,6,1,6);
  double **derivs = dmatrix(1,6,1,6);
  /* Transform the orbit basis and get the deriv. matrix */
  pbasis_to_bary(&p, &xv, derivs);

  /* Map the covariance matrix to new basis */
  covar_map(covar, derivs, covar_xyz,6,6);

  /* Get partial derivative matrix from xyz to aei */
  aei_derivs(&xv, derivs);
  /* Map the covariance matrix to new basis */
  covar_map(covar_xyz, derivs, covar_aei,6,6);

  /* Transform xyz basis to orbital parameters */
  orbitElements(&xv, &orbit);

    /* Print out the results, with comments */
  printf("# Barycentric osculating elements in ICRS at epoch %.1f:\n",jd0);
  printf("#    a            e       i      Node   Arg of Peri   Time of Peri\n");
  printf("%12.6f  %9.6f  %8.3f %8.3f  %8.3f %11.3f\n",
	 orbit.a, orbit.e, orbit.i, orbit.lan, orbit.aop, orbit.T);
  printf("+-%10.6f  %9.6f  %8.3f %8.3f  %8.3f %11.3f\n",
	 sqrt(covar_aei[1][1]),
	 sqrt(covar_aei[2][2]),
	 sqrt(covar_aei[3][3])/DTOR,
	 sqrt(covar_aei[4][4])/DTOR,
	 sqrt(covar_aei[5][5])/DTOR,
	 sqrt(covar_aei[6][6])/DAY);
  printf("# aei covariance matrix:\n");

/* Hack to return the covariance matrix in a format that the Python wrapper can understand */
   double cov1d[36];
   flatten_cov(covar, 6, cov1d);

  print_matrix(stdout,covar_aei,6,6);

  free_dmatrix(covar,1,6,1,6);
  free_dmatrix(covar_xyz,1,6,1,6);
  free_dmatrix(covar_aei,1,6,1,6);
} 

