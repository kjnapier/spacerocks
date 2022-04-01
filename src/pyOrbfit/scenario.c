/* 	$Id: scenario.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: scenario.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* scenario.c - Read a file containing a/b/g orbit and 
 * initial dates of observation.  Make mock observations and
 * re-fit.  See how posn uncertainty grows; make new mock observations
 * when the posn uncertainty grows to a certain level.
 * 6/14/00 gmb
 */
#include "orbfit.h"

#define TIMESTEP 1.1

char   *help[] = {
    "scenario: Read a file containing a/b/g orbit fit, and input",
    "       dates of initial observations.  Then watch uncertainties",
    "       grow and `reobserve' when they are big.",
    " usage: scenario <abgfile> <threshold> <timespan>",
    " abgfile    is name of file with orbit info (from fit_radec)",
    " threshold  is maximal allowed position uncertainty (arcsec)",
    " timespan   is how long (years) to run the simulation",
    " input      is list of initial observations, one per line.",
    "            Note observatory and errors will be taken from",
    "            the last input line."
    " output     is log of uncertainties & observations",
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

void
print_radec(OBSERVATION *obs, FILE *fptr)
{
  /* Print out new observation */      
  /* Now transform to RA/DEC, via ecliptic*/
  double ra, dec, lat, lon;
  char rastring[20],decstring[20];
  void deghms(double degr,
	      char *outbuff);
  void degdms(double degr,
	      char *outbuff); 
  proj_to_ec(obs->thetax,obs->thetay,
	     &lat, &lon,
	     lat0, lon0, NULL);
  ec_to_eq(lat, lon, &ra, &dec, NULL);
  ra /= DTOR;
  if (ra<0.) ra+= 360.;
  dec /= DTOR;
  deghms(ra,rastring);
  degdms(dec,decstring);
  fprintf(fptr,"%.4f %s %s %5.2f %d\n", obs->obstime/DAY + jd0,
	 rastring,decstring,obs->dthetax/ARCSEC,
	 obs->obscode);
  return;
}


int
main(int argc, char *argv[])
{
  PBASIS ptrue, pfit;
  OBSERVATION	futobs, obsarray[MAXOBS], *obs;
  int nobs;
  double tt;

  ORBIT  orbit;
  XVBASIS xv;

  double **covar, **covar_aei, **covar_xyz;
  double chisq;
  int dof;

  char	inbuff[256];
  double **sigxy,a,b,PA,**derivs;
  double lat,lon,**covecl;
  double ra,dec, **coveq;
  double xx,yy,xy,bovasqrd,det;
  double now, threshold, timestart, timespan;
  int i, refit;

  if (argc!=4 || *argv[1]=='^') print_help();


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
  covar_xyz = dmatrix(1,6,1,6);
  covar_aei = dmatrix(1,6,1,6);
  covecl = dmatrix(1,2,1,2);
  coveq = dmatrix(1,2,1,2);

  if (read_abg(argv[1],&ptrue,covar)) {
    fprintf(stderr, "Error input alpha/beta/gamma file %s\n",argv[1]);
    exit(1);
  }

  if ((threshold = atof(argv[2]))<=0.) {
    fprintf(stderr,"Bad observation threshold %s\n",argv[2]);
    print_help();
  }

  if ((timespan = atof(argv[3]))<=0.) {
    fprintf(stderr,"Bad time span %s\n",argv[3]);
    print_help();
  }

  /* Read in the initial observation times */
  nobs = 0;
  while (fgets_nocomment(inbuff, 255, stdin, stdout)!=NULL) {
    obs = &(obsarray[nobs]);
    if (scan_observation(inbuff, obs)) exit(1);
    /* Set time to years after jd0, don't care about coordinates here*/
    obs->obstime = (obs->obstime-jd0)*DAY;
    /* Calculate the position of Earth at this time to avoid doing
     * it many times later: */
    earth3d(obs->obstime, obs->obscode,
	    &(obs->xe),&(obs->ye),&(obs->ze));
    /* replace with faked position */
    fake_observation(&ptrue, obs);
    print_radec(obs,stdout);
    nobs++;
  }

  /*Save site and error limits to use in future */
  futobs.obscode = obsarray[nobs-1].obscode;
  futobs.dthetax = obsarray[nobs-1].dthetax;
  futobs.dthetay = obsarray[nobs-1].dthetay;

  timestart = obsarray[0].obstime;

  refit = 1;
  for (now = obsarray[nobs-1].obstime; now < timespan; 
       now = timestart+TIMESTEP*(now-timestart)) {

    /*first time through or with large uncertainty, refit orbit*/
    if (refit) {
      refit = 0;
      fit_observations(obsarray, nobs, &pfit, covar, &chisq, &dof);
      /*Spit out a, e, and their scaled uncertainties*/
      pbasis_to_bary(&pfit, &xv, derivs);
      covar_map(covar, derivs, covar_xyz,6,6);
      orbitElements(&xv, &orbit);
      aei_derivs(&xv, derivs);
      covar_map(covar_xyz, derivs, covar_aei,6,6);

      /* Compute a, b, theta of error ellipse in a-e plane when scaled to a & e*/
      xx = covar_aei[1][1]/(orbit.a*orbit.a);
      xy = covar_aei[1][2]/(orbit.a*orbit.e);
      yy = covar_aei[2][2]/(orbit.e*orbit.e);
      PA = 0.5 * atan2(2.*xy,(xx-yy)) * 180./PI;	/*go right to degrees*/
      bovasqrd  = (xx+yy-sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) 
	/ (xx+yy+sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) ;
      det = xx*yy-xy*xy;
      b = pow(det*bovasqrd,0.25);
      a = pow(det/bovasqrd,0.25);

      printf("***New orbit fit: a,e: %.4f %.4f %.4f %.4f %.2f\n",
	     orbit.a, orbit.e, a, b, PA);
    }

    /* Predict position & examine uncertainty */
    futobs.obstime = now;
    futobs.xe = -999.;		/* Force evaluation of earth3d */

    predict_posn(&pfit,covar,&futobs,sigxy);
    /* Compute a, b, theta of error ellipse for output */
    xx = sigxy[1][1];
    yy = sigxy[2][2];
    xy = sigxy[1][2];
    PA = 0.5 * atan2(2.*xy,(xx-yy)) * 180./PI;	/*go right to degrees*/
    /* Adjust for PA to be N through E, */
    PA = PA-90;
    if (PA<-90.) PA += 180.;
    bovasqrd  = (xx+yy-sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) 
      / (xx+yy+sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) ;
    det = xx*yy-xy*xy;
    b = pow(det*bovasqrd,0.25);
    a = pow(det/bovasqrd,0.25);

    printf("%7.3f %9.2f  %9.2f   %8.2f %8.2f %7.2f\n",
	   now, futobs.thetax/ARCSEC, futobs.thetay/ARCSEC,
	   a/ARCSEC,b/ARCSEC,PA); 

    /* If the uncertainty is above threshold, generate 
     * a new pair of observations*/
    if (a/ARCSEC>threshold) {
      obsarray[nobs].obstime = now;
      obsarray[nobs].obscode = futobs.obscode;
      obsarray[nobs].dthetax = futobs.dthetax;
      obsarray[nobs].dthetay = futobs.dthetay;
      obsarray[nobs].xe = -999.;
      fake_observation(&ptrue, &(obsarray[nobs]));
      print_radec(obs,stdout);
      nobs++;

      obsarray[nobs].obstime = now+DAY;
      obsarray[nobs].obscode = futobs.obscode;
      obsarray[nobs].dthetax = futobs.dthetax;
      obsarray[nobs].dthetay = futobs.dthetay;
      obsarray[nobs].xe = -999.;
      fake_observation(&ptrue, &(obsarray[nobs]));
      print_radec(obs,stdout);
      nobs++;

      refit = 1;
    }
  }
   
  free_dmatrix(covar,1,6,1,6);
  free_dmatrix(covar_xyz,1,6,1,6);
  free_dmatrix(covar_aei,1,6,1,6);
  free_dmatrix(derivs,1,6,1,6);
  exit(0);
}
