/* 	$Id: scenario2.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: scenario2.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* scenario2.c - Read a file containing a/b/g orbit and 
 * initial dates of observation.  Make mock observations and
 * re-fit.  
 * This time, check mock observations for elongation constraint.
 * Also choose the possible observation which minimizes sigma_a.
 * 6/16/00 gmb
 */
#include "orbfit.h"

/* maximum elongation at which observation is possible*/
#define MIN_ELONG (90.*DTOR)
#define MONTH (29.*DAY)

char   *help[] = {
    "scenario2: Read a file containing a/b/g orbit fit, and input",
    "       dates of initial observations.  Then watch uncertainties",
    "       grow and `reobserve' when sigma_a is minimized.",
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
  PBASIS ptrue, pfit, ptry;
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
  double now, threshold, timestart, timespan, timeend;
  double besttime, bestvara, elong;
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
    if (nobs==0) timestart = obs->obstime;
    /* Calculate the position of Earth at this time to avoid doing
     * it many times later: */
    earth3d(obs->obstime, obs->obscode,
	    &(obs->xe),&(obs->ye),&(obs->ze));
    /* replace with faked position */
    fake_observation(&ptrue, obs);
    print_radec(obs,stdout);
    printf("Reobserve time %.4f opp. angle %.2f\n", obs->obstime-timestart,
	   opposition_angle(obs)/DTOR);
    nobs++;
  }

  /* Get an orbit fit*/
  fit_observations(obsarray, nobs, &pfit, covar, &chisq, &dof, NULL);

  /*Save site and error limits to use in future */
  obsarray[nobs].obscode = obsarray[nobs-1].obscode;
  obsarray[nobs].dthetax = obsarray[nobs-1].dthetax;
  obsarray[nobs].dthetay = obsarray[nobs-1].dthetay;

  timeend = timestart + timespan;
  bestvara = 1e10;
  besttime = 0.;

  for (now = obsarray[nobs-1].obstime; now <= timeend;  ) {

    /*Advance observing time */
    if ( now-timestart < 2*MONTH) now+= DAY;
    else now+=MONTH;

    /*Get position and x uncertainty*/
    /* Predict position & examine uncertainty */
    obsarray[nobs].obstime = now;
    obsarray[nobs].xe = -999.;		/* Force evaluation of earth3d */

    predict_posn(&pfit,covar,&obsarray[nobs],sigxy);


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

    /* See what uncertainty in a results from this obs.*/
    fit_observations(obsarray, nobs+1, &ptry, covar_aei, &chisq, &dof, NULL);

    pbasis_to_bary(&ptry, &xv, derivs);
    covar_map(covar_aei, derivs, covar_xyz,6,6);
    orbitElements(&xv, &orbit);
    aei_derivs(&xv, derivs);
    covar_map(covar_xyz, derivs, covar_aei,6,6);
    
    /* Report results of this test to output */
    printf("%7.4f %.2f %.2f %.2f %9.6f %.3g %5.1f  %5.1f %d\n", now-timestart,
	   obsarray[nobs].thetax/ARCSEC, obsarray[nobs].thetay/ARCSEC,
	   a/ARCSEC, sqrt(covar_aei[1][1]), sqrt(covar_aei[2][2]),
	   elong/DTOR, opposition_angle(&obsarray[nobs])/DTOR, nobs);
    
    /* Go to next date if Sun is too close */
    elong = elongation(&obsarray[nobs]);
    if (elong < MIN_ELONG) continue;

    if (a > threshold*ARCSEC) {
      /* object is lost.  See what prior observation was best
       * for constraining a, then add it and refit orbit, continue
       */
      if (besttime<=0) {
	fprintf(stderr,"Object is lost before observation\n");
	exit(1);
      }
      obsarray[nobs].obstime = besttime;
      obsarray[nobs].xe = -999.;
      fake_observation(&ptrue, &obsarray[nobs]);
      printf("Reobserve time %.4f opp. angle %.2f\n", besttime-timestart,
	    opposition_angle(&obsarray[nobs])/DTOR);
      print_radec(&obsarray[nobs],stdout);
      nobs++;
      fit_observations(obsarray, nobs, &pfit, covar, &chisq, &dof, NULL);
      now = besttime;
      besttime = 0.;
      bestvara = 1e10;
      obsarray[nobs].obscode = obsarray[nobs-1].obscode;
      obsarray[nobs].dthetax = obsarray[nobs-1].dthetax;
      obsarray[nobs].dthetay = obsarray[nobs-1].dthetay;
    } else {

      /* Is this the best observation so far? */
      if (covar_aei[1][1] < bestvara) {
	bestvara = covar_aei[1][1];
	besttime = now;
      }
    }
  }

  /*Do a final observation*/
  obsarray[nobs].obstime = besttime;
  obsarray[nobs].xe = -999.;
  fake_observation(&ptrue, &obsarray[nobs]);
  printf("Reobserve time %.4f opp. angle %.2f\n", besttime-timestart,
	 opposition_angle(&obsarray[nobs])/DTOR);
  print_radec(&obsarray[nobs],stdout);
  nobs++;
  fit_observations(obsarray, nobs, &pfit, covar_aei, &chisq, &dof, NULL);

  pbasis_to_bary(&pfit, &xv, derivs);
  covar_map(covar_aei, derivs, covar_xyz,6,6);
  orbitElements(&xv, &orbit);
  aei_derivs(&xv, derivs);
  covar_map(covar_xyz, derivs, covar_aei,6,6);
  printf("Final uncertainty in a: %.3g\n",
	 sqrt(covar_aei[1][1]));

   
  free_dmatrix(covar,1,6,1,6);
  free_dmatrix(covar_xyz,1,6,1,6);
  free_dmatrix(covar_aei,1,6,1,6);
  free_dmatrix(derivs,1,6,1,6);
  exit(0);
}
