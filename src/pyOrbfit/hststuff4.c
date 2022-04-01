/* 	$Id: hststuff4.c,v 2.0 2001/09/14 19:05:06 garyb Exp $ */

#ifndef lint
static char vcid[] = "$Id: hststuff4.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */

#include "orbfit.h"
OBSERVATION obsarray[MAXOBS];
int	nobs;

/*suffix for orbit-fit files*/
#define SUFFIX ".abg"
#define OBSCODE_HST 2000
#define TIME_STEP   0.2	  /*time step (days) to locate oppangle*/
#define MAX_EPOCHS  100
#define LAT_MAX     5     /*only objects that are this close to ecliptic*/
#define PHASE_STEPS 20	   /*number of phases for each object's orbit*/

#define ZA_LIMIT (90.*DTOR)	/*maximum zenith angle for observations*/
#define NORBITS  10	/*number of orbits to observe*/
#define NFIELDS  4      /*number of fields to observe in cycle*/
#define EXP_PER_ORBIT 6	/*exposures per orbit*/
#define ASTROM_ERROR  0.004	/*uncertainty per visit */
#define VIS_TIME_STEP (1./24./60.) /*1-minute time steps to check visibility*/

#define OBSCODE_GRND 807	/*Where the followup observation will be*/
#define ERROR_GRND  0.2	   /*Astrometric error on followup*/
#define OPP_GRND    225.   /*opp angle to do pair of ground observations*/

char   *help[] = {
  "hststuff4:  Show orbit uncertainties based upon a series of",
  "            HST observations, plus perhaps an additional ground-based",
  "            observation later on.",
  "            Try each object at several points in its orbit."
  " usage: hststuff4 oppangle day1 day2 ... dayn",
  " oppangle   is opposition angle for first observation, which will",
  "            serve as zeropoint for positions.",
  " dayn       are days after first observation at which to calculate",
  " input      is list of objects.  It is assumed that a file",
  "            with the name <object>.abg is present in current directory",
  " output     are predicted RA & Dec plus error ellipse & elongation.",
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

/* routine to advance jd until oppangle is reached */
void
seek_oppangle(double oppangle,
	      double *jd,
	      PBASIS *p,
	      int    obscode)
{
  OBSERVATION futobs;
  double nextopp, difflast, diffnext, lastopp;
  
  futobs.obscode = obscode;
  futobs.obstime = (*jd-jd0)*DAY;
  futobs.xe = -999.;		/* Force evaluation of earth3d */
  predict_posn(p,NULL,&futobs,NULL);
  lastopp = opposition_angle(&futobs)/DTOR;

  do {
    *jd += TIME_STEP;
    futobs.obstime = (*jd-jd0)*DAY;
    futobs.xe = -999.;		/* Force evaluation of earth3d */
    predict_posn(p,NULL,&futobs,NULL);
    nextopp = opposition_angle(&futobs)/DTOR;
    
    difflast = oppangle - lastopp;
    while (difflast < -180.) difflast+=360.;
    while (difflast > 180.) difflast-=360.;
    diffnext = oppangle - nextopp;
    while (diffnext < -180.) diffnext+=360.;
    while (diffnext > 180.) diffnext-=360.;

    lastopp = nextopp;
  } while ( diffnext != 0. && 
	    (diffnext*difflast>0. || fabs(diffnext)>120.));

  return;
}

int
main(int argc, char *argv[])
{
  PBASIS p, pfit;
  ORBIT o, ofit;
  XVBASIS xv;

  OBSERVATION	futobs;
  OBSERVATION *oo;
  FILE *abgfile;
  char	inbuff[256],objname[256],fname[256];

  double **covar,**sigxy;
  double lat,lon;
  int i,nfields;
  double x[3], v[3], range;

  double oppangle, epochs[MAX_EPOCHS], observe_day, period;
  double x0, y0;
  int nepochs, iphase, iexpose, iorbit;

  int dof,ii;
  double chisq, dbary, ddbary, vis_start, vis_end, jd;
  double  **covar_abg, **covar_xyz, **derivs, **covar_aei;
  double launch;	/*Date to start hunting for oppangle*/

  if (argc<2) print_help();
  oppangle = atof(argv[1]);
  epochs[0]=0.;
  for (i=2; i<argc; i++) epochs[i-1] = atof(argv[i]);
  nepochs = argc-1;

  futobs.obscode = OBSCODE_HST;

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

  printf(
"#   ID  phi  distance            a                   e                     i   "
"    posn        d             a              e              i \n"
"#                        real   fit     err    real   fit     err     fit  err  NEW OBSERVATION:\n");
 
  covar_abg = dmatrix(1,6,1,6);
  covar_xyz = dmatrix(1,6,1,6);
  covar_aei = dmatrix(1,6,1,6);
  derivs = dmatrix(1,6,1,6);
  sigxy = dmatrix(1,2,1,2);
  covar = dmatrix(1,6,1,6);

  /* Loop through input containing names of objects */
  while ( fgets_nocomment(inbuff,255,stdin,NULL)!=NULL) {
    sscanf(inbuff," %s ",objname);
    strcpy(fname,objname);
    strcat(fname,SUFFIX);

    if (read_abg(fname,&p,covar)) {
      fprintf(stderr, "Error reading alpha/beta/gamma file %s\n",fname);
      exit(1);
    }

    /* look at several points around the orbital phase */
    /* get orbital elements for this one */
    pbasis_to_bary(&p, &xv, NULL);
    orbitElements(&xv, &o);

    period = pow(o.a, 1.5);
    launch = jd0;

    for (iphase=0; iphase<PHASE_STEPS; iphase++) {
      /* Advance the object through its orbit by moving back time
       * of perihelion (which is in days) */
      o.T -= period/ (PHASE_STEPS*DAY);
      /* Set up new orbit system based on this */
      elements_to_pbasis(&o, launch, futobs.obscode, &p);
      
      /* Now look for the day at which opp. angle meets spec */
      observe_day = launch;
      seek_oppangle(oppangle, &observe_day, &p, futobs.obscode);

      futobs.obstime = (observe_day-jd0)*DAY;
      futobs.xe = -999.;		/* Force evaluation of earth3d */
      predict_posn(&p,NULL,&futobs,NULL);

      /* Now transform ecliptic*/
      proj_to_ec(futobs.thetax,futobs.thetay,
		 &lat, &lon,
		 lat0, lon0, NULL);

      /* Only keep objects close to ecliptic */
      if (fabs(lat) > LAT_MAX*DTOR) continue;

      /* Now reset the orbit and coordinate system to start here */
      elements_to_pbasis(&o, observe_day, futobs.obscode, &p);
      nobs = 0;

      /* accumulate observations at all epochs */
      for (i=0; i<nepochs; i++) {
	/* Advance time until object is NOT visible */
	vis_end = observe_day+epochs[i];
	do {
	  vis_end+= VIS_TIME_STEP;
	  futobs.obstime = (vis_end-jd0)*DAY;
	  futobs.xe = -999.;		/* Force evaluation of earth3d */
	  predict_posn(&p,NULL,&futobs,NULL);
	} while (zenith_angle(&futobs) <= ZA_LIMIT) ;

	for (iorbit=0; iorbit<NORBITS*NFIELDS; iorbit++) {
	  /* Then advance to first time it IS visible. */
	  vis_start = vis_end;
	  do {
	    vis_start+= VIS_TIME_STEP;
	    futobs.obstime = (vis_start-jd0)*DAY;
	    futobs.xe = -999.;		/* Force evaluation of earth3d */
	    predict_posn(&p,NULL,&futobs,NULL);
	  } while (zenith_angle(&futobs) > ZA_LIMIT) ;
	  
	  /* Now find end of visibility period */
	  vis_end = vis_start;
	  do {
	    vis_end+= VIS_TIME_STEP;
	    futobs.obstime = (vis_end-jd0)*DAY;
	    futobs.xe = -999.;		/* Force evaluation of earth3d */
	    predict_posn(&p,NULL,&futobs,NULL);
	  } while (zenith_angle(&futobs) <= ZA_LIMIT) ;
	  vis_end -= VIS_TIME_STEP;
	  
	  /* Only proceed if we are observing this field on this orbit*/
	  if (iorbit%NFIELDS !=0) continue;
	  
	  /* Find middle of each exposure time */
	  jd = vis_start + (vis_end - vis_start)/(2.*EXP_PER_ORBIT);
	  for (iexpose=0; iexpose<EXP_PER_ORBIT; iexpose++, 
		 jd+=(vis_end - vis_start)/EXP_PER_ORBIT) {
	    oo = &obsarray[nobs];
	    oo->obstime = (jd-jd0)*DAY;
	    oo->xe = -999.;		/* Force evaluation of earth3d */
	    oo->obscode = OBSCODE_HST;
	    predict_posn(&p,NULL,oo,NULL);
	    oo->dthetax = oo->dthetay =  
	      ASTROM_ERROR*sqrt(NORBITS*EXP_PER_ORBIT)*ARCSEC;
	    nobs++;

	    /*if (iphase==0) {
	      double ra, dec;
	      char rastring[40], decstring[40];
	      proj_to_ec(oo->thetax,oo->thetay,
			 &lat, &lon,
			 lat0, lon0, NULL);
	      ec_to_eq(lat, lon, &ra, &dec, NULL);
	      deghms(ra/DTOR,rastring);
	      degdms(dec/DTOR,decstring);
	      fprintf(stderr,"%f %s %s %.3f %d\n",jd, rastring, decstring,
		      oo->dthetax/ARCSEC,
		      OBSCODE_HST);
	    }*/

	  } /*exposure loop*/
	} /*orbit loop*/
      } /*observation-epoch loop */

      /* Now fit the observations and report results */
      ii=fit_observations(obsarray, nobs, &pfit, covar, &chisq, &dof, NULL);
      /*fprintf(stderr,"#     lat0       lon0       xBary     yBary      zBary   JD0\n");
	fprintf(stderr,"%12.7f %12.7f %10.7f %10.7f %10.7f  %.6f\n",
	lat0/DTOR,lon0/DTOR,xBary,yBary,zBary,jd0);
	fprintf(stderr,"# Chi-squared of fit: %.2f DOF: %d\n",chisq,dof);
	fprintf(stderr,"%11.8f %11.8f %11.8f %11.8f %11.8f %11.8f\n",
	pfit.a,pfit.adot,pfit.b,pfit.bdot, pfit.g, pfit.gdot);
	print_matrix(stderr,covar,6,6);
      */
      pbasis_to_bary(&pfit, &xv, derivs);
      orbitElements(&xv, &ofit);
      /* Note that below is distance at jd0=abg file zpoint, NOT at
       * first observation, which is some part of a year later.
       */
      dbary = sqrt(xBary*xBary + yBary*yBary + pow(zBary-1/pfit.g,2.));
      ddbary = dbary*dbary*sqrt(covar[5][5]);

      /* Get elements and covariance matrix thereof */
      covar_map(covar, derivs, covar_xyz,6,6);
      aei_derivs(&xv, derivs);
      covar_map(covar_xyz, derivs, covar_aei,6,6);
    
      /*Print out some info & uncertainties */
      printf("%8s %2d %5.2f +- %4.2f %1d %5.2f %5.2f +- %4.2f %5.3f %5.3f +- %5.3f"
	     " %7.1f +- %4.1f",
	     objname, iphase, 
	     dbary, ddbary, ii,
	     o.a, ofit.a, sqrt(covar_aei[1][1]),
	     o.e, ofit.e, sqrt(covar_aei[2][2]),
	     ofit.i, sqrt(covar_aei[3][3])/DTOR);

      /* Now get the positional uncertainty at a future date, and
       * refit after two extra measurements from ground
       */
      oo = &obsarray[nobs];
      nobs++;
      jd = observe_day + OPP_GRND - oppangle;
      oo->obstime = (jd-jd0)*DAY;
      oo->obscode = OBSCODE_GRND;
      oo->xe = -999.;
      predict_posn(&pfit,covar,oo,sigxy);
      /*if (iphase==0) {
	fprintf(stderr,"Prediction on %.2f, posn %.2f %.2f unc %.2f %.2f\n",
		    jd, oo->thetax/ARCSEC, oo->thetay/ARCSEC,
		    sqrt(sigxy[1][1])/ARCSEC, sqrt(sigxy[2][2])/ARCSEC);
	print_matrix(stderr,covar,6,6);
      }*/

      { 
	/* Compute a, of error ellipse for output */
	double xx, yy, xy, bovasqrd, det, a;
	xx = sigxy[1][1];
	yy = sigxy[2][2];
	xy = sigxy[1][2];
	bovasqrd  = (xx+yy-sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) 
	  / (xx+yy+sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) ;
	det = xx*yy-xy*xy;
	a = pow(det/bovasqrd,0.25);
	printf(" %5.1f", a/ARCSEC);
      }


      /* Now re-fit with two additional observations */
      predict_posn(&p,NULL,oo,NULL);
      oo->dthetax = oo->dthetay = ERROR_GRND * ARCSEC;
      oo = &obsarray[nobs];
      nobs++;
      oo->obstime = obsarray[nobs-2].obstime + DAY;
      oo->obscode = OBSCODE_GRND;
      oo->xe = -999.;
      predict_posn(&p,NULL,oo,NULL);
      oo->dthetax = oo->dthetay = ERROR_GRND * ARCSEC;

      ii = fit_observations(obsarray, nobs, &pfit, covar, &chisq, &dof, NULL);
      pbasis_to_bary(&pfit, &xv, derivs);
      orbitElements(&xv, &ofit);

      dbary = sqrt(xBary*xBary + yBary*yBary + pow(zBary-1/pfit.g,2.));
      ddbary = dbary*dbary*sqrt(covar[5][5]);

      /* Get elements and covariance matrix thereof */
      covar_map(covar, derivs, covar_xyz,6,6);
      aei_derivs(&xv, derivs);
      covar_map(covar_xyz, derivs, covar_aei,6,6);
    
      /*Print out some info & uncertainties */
      printf(" %5.2f +- %4.2f %1d %5.2f +- %4.2f %5.3f +- %5.3f"
	     " %7.1f +- %4.1f \n",
	     dbary, ddbary, ii,
	     ofit.a, sqrt(covar_aei[1][1]),
	     ofit.e, sqrt(covar_aei[2][2]),
	     ofit.i, sqrt(covar_aei[3][3])/DTOR);

    } /*phase loop */
  } /*object loop*/
  exit(0);
}
