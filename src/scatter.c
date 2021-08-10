/* 	$Id: scatter.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: scatter.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* scatter.c - Read a file containing a/b/g orbit fit, and
 * a list of observations.  Fit orbit, make prediction & error,
 * then randomize the input observations and see what the distrib
 * of predicted values is w.r.t. true value.
 * 6/14/00 gmb
 */
#include "orbfit.h"

char   *help[] = {
    "scatter: Read a true orbit fit, and list of observations, then",
    "       generate random data to check the scatter",
    " usage: predict [-m mpc_error] [-j JPL_file] [-o observatory_file] [-v]",
    "                <abgfile> <ntrials> <xy_or_ae>",
    "  mpc_error  is arcsec uncertainty to apply to MPC-format",
    "             observations.  Default is 0.2",
    "  JPL_file   is binary ephemeris file.  Default is binEphem.423, or",
    "             a file specified by environment variable ORBIT_EPHEMERIS",
    "  observatory_file   is file with observatory site data.  Default is",
    "             observatories.dat, or a file specified by environment",
    "             variable ORBIT_OBSERVATORIES",
    " abgfile    is name of file with orbit info (from fit_radec)",
    " ntrials    is number of random trials to do",
    " xy_or_ae   is 0 to predict x&y locations (at last obs), 1 for a/e",
    " input      are observations in MPC or my format.",
    " output     residuals plus chi-squared values",
    0
};

void
print_help(void)
{
  int	i;
  for (i = 0; help[i] != 0; i++)
    fprintf (stderr, "%s\n", help[i]);
  exit (1);
}

extern void 
deghms(double degr,
       char *outbuff);
extern void 
degdms(double degr,
       char *outbuff);

int
main(int argc, char *argv[])
{
  PBASIS ptrue, pfit;
  OBSERVATION	obsarray[MAXOBS], testobs, *obs;
  XVBASIS xv;
  ORBIT orbit;
  int nobs;
  struct date_time dt;
  char	inbuff[256],rastring[20],decstring[20];
  double **covar,**sigxy,a,b,PA,**derivs, **covar_aei, **covar_xyz;
  double lat,lon,**covecl;
  double ra,dec, **coveq;
  double dx, dy, chisq, elat, elon;
  double xx,yy,xy,bovasqrd,det;
  int i,nfields, dof;

  int doing_ae, ntrials, itrial;
  double truex, truey, sumx, sumy, sumxx, sumyy, sumxy, dev;

  int iarg=1;
  if (argc>1 && *argv[1]=='^') print_help();
  if (read_options(&iarg, argc, argv)) print_help();
  if (argc-iarg!=3) print_help();


  /* echo the command line to output */
  printf("#");
  for (i=0; i<argc; i++) printf(" %s",argv[i]);
  {
#include <time.h>
    time_t timettt;
    time(&timettt);
    /* note that ctime returns string with newline at end */
    printf("\n#---%s",ctime(&timettt));
  }

  sigxy = dmatrix(1,2,1,2);
  derivs = dmatrix(1,6,1,6);
  covar = dmatrix(1,6,1,6);
  covecl = dmatrix(1,2,1,2);
  coveq = dmatrix(1,2,1,2);
  covar_xyz = dmatrix(1,6,1,6);
  covar_aei = dmatrix(1,6,1,6);

  if (read_abg(argv[iarg],&ptrue,covar)) {
    fprintf(stderr, "Error input alpha/beta/gamma file %s\n",argv[iarg]);
    exit(1);
  }
  iarg++;

  ntrials = atoi(argv[++iarg]);
  doing_ae = atoi(argv[iarg]);
  if (ntrials<=0 || doing_ae<0 || doing_ae>1) {
    fprintf(stderr,"Error specifying ntrials or xy/ae\n");
    print_help();
  }
  /* Read all observation conditions from stdin, fill with correct posns*/
  nobs = 0;
  while (fgets_nocomment(inbuff, 255, stdin, stdout)!=NULL) {
    obs = &(obsarray[nobs]);
    if (scan_observation(inbuff, obs)) exit(1);
    /* Set time to years after jd0, rotate to tangent plane coords */
    obs->obstime = (obs->obstime-jd0)*DAY;
    /* Calculate the position of Earth at this time to avoid doing
     * it many times later: */
    earth3d(obs->obstime, obs->obscode,
	    &(obs->xe),&(obs->ye),&(obs->ze));
    kbo2d(&ptrue,obs,&(obs->thetax),NULL,&(obs->thetay),NULL);
    nobs++;
  }

  /* If we are going to measure xy scatter, the last observation was
   * the "test" time.  Otherwise we can get the true a & e from the abg file.
   */
  if (doing_ae) {
    pbasis_to_bary(&ptrue, &xv, NULL);
    orbitElements(&xv, &orbit);
    truex = orbit.a;
    truey = orbit.e;
  } else {
    nobs--;
    testobs.obstime = obsarray[nobs].obstime;
    testobs.obscode = obsarray[nobs].obscode;
    testobs.xe = obsarray[nobs].xe;
    testobs.ye = obsarray[nobs].ye;
    testobs.ze = obsarray[nobs].ze;
    truex = obsarray[nobs].thetax;
    truey = obsarray[nobs].thetay;
  }

  /* Fit the exact data to get its errors */
  fit_observations(obsarray, nobs, &pfit, covar, &chisq, &dof, NULL);
  
  /*Get and report errors appropriate to this case: */
  if (doing_ae) {
    /* Transform the orbit basis and get the deriv. matrix */
    pbasis_to_bary(&pfit, &xv, derivs);
    covar_map(covar, derivs, covar_xyz,6,6);
    aei_derivs(&xv, derivs);
    covar_map(covar_xyz, derivs, covar_aei,6,6);

    /* Transform xyz basis to orbital parameters */
    orbitElements(&xv, &orbit);
    printf("# Est  a, e, variance components:\n");
    printf("%10.6f %10.6f %.5g %.5g %.5g\n",truex, truey,
	   covar_aei[1][1], covar_aei[2][2], covar_aei[1][2]);
    xx = covar_aei[1][1];
    yy = covar_aei[2][2];
    xy = covar_aei[1][2];
  } else {
    predict_posn(&pfit,covar,&testobs,sigxy);
    printf("# Est  x, y, variance components:\n");
    printf("%10.6f %10.6f %.5g %.5g %.5g\n",truex/ARCSEC, truey/ARCSEC,
	   sigxy[1][1]/(ARCSEC*ARCSEC), sigxy[2][2]/(ARCSEC*ARCSEC), 
	   sigxy[1][2]/(ARCSEC*ARCSEC));
    xx = sigxy[1][1]/(ARCSEC*ARCSEC);
    yy = sigxy[2][2]/(ARCSEC*ARCSEC);
    xy = sigxy[1][2]/(ARCSEC*ARCSEC);

  }
  /* Compute a, b, theta of error ellipse for output */
  PA = 0.5 * atan2(2.*xy,(xx-yy)) * 180./PI;	/*go right to degrees*/
  if (PA<-90.) PA += 180.;
  bovasqrd  = (xx+yy-sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) 
    / (xx+yy+sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) ;
  det = xx*yy-xy*xy;
  b = pow(det*bovasqrd,0.25);
  a = pow(det/bovasqrd,0.25);
  printf("Error ellipse a,b, PA (x-axis=0): %g %g %.2f\n",a,b,PA);

  /* Now start the random trials. */
  sumx = sumy = sumxx = sumyy = sumxy = 0.;
  for (itrial = 0; itrial<ntrials; itrial++) {
    /*Generate & fit randomized observations */
    for (i=0; i<nobs; i++) {
          fake_observation(&ptrue, &(obsarray[i]));
    }
    fit_observations(obsarray, nobs, &pfit, covar, &chisq, &dof, NULL);

    /* Get the desired quantities from this realization:*/
    if (doing_ae) {
      /* Transform the orbit basis and get the deriv. matrix */
      pbasis_to_bary(&pfit, &xv, derivs);
      covar_map(covar, derivs, covar_xyz,6,6);
      aei_derivs(&xv, derivs);
      covar_map(covar_xyz, derivs, covar_aei,6,6);

      /* Transform xyz basis to orbital parameters */
      orbitElements(&xv, &orbit);
      dx = orbit.a - truex;
      dy = orbit.e - truey;
      dev = dx*dx*covar_aei[2][2] + dy*dy*covar_aei[1][1]
	- 2*dx*dy*covar_aei[1][2];
      dev /= covar_aei[1][1]*covar_aei[2][2]-covar_aei[1][2]*covar_aei[2][1];
      printf("Fit result: %5.2f %9.6f %9.6f %5.2f\n", chisq, dx, dy, dev);
    } else {
      predict_posn(&pfit,covar,&testobs,sigxy);
      dx = testobs.thetax - truex;
      dy = testobs.thetay - truey;
      dev = dx*dx*sigxy[2][2] + dy*dy*sigxy[1][1]
	- 2*dx*dy*sigxy[1][2];
      dev /= sigxy[1][1]*sigxy[2][2]-sigxy[1][2]*sigxy[2][1];
      dx /= ARCSEC;
      dy /= ARCSEC;
      printf("Fit result: %5.2f %9.3f %9.3f %5.2f\n", chisq, 
	     dx, dy, dev);
    }

    /* accumulate errors */
    sumx += dx;
    sumy += dy;
    sumxx += dx*dx;
    sumyy += dy*dy;
    sumxy += dx*dy;
    
  }

  sumx /= ntrials;
  sumy /= ntrials;
  sumxx /= ntrials;
  sumyy /= ntrials;
  sumxy /= ntrials;
  xx = sumxx - sumx*sumx;
  xy = sumxy - sumx*sumy;
  yy = sumyy - sumy*sumy;
  printf("Actual mean, variance: %g %g %g %g %g\n",
	 sumx, sumy, xx, yy, xy);
  PA = 0.5 * atan2(2.*xy,(xx-yy)) * 180./PI;	/*go right to degrees*/
  if (PA<-90.) PA += 180.;
  bovasqrd  = (xx+yy-sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) 
    / (xx+yy+sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) ;
  det = xx*yy-xy*xy;
  b = pow(det*bovasqrd,0.25);
  a = pow(det/bovasqrd,0.25);
  printf("Error ellipse a,b, PA (x-axis=0): %g %g %.2f\n",a,b,PA);

  exit(0);
}
