/* 	$Id: hststuff3.c,v 2.0 2001/09/14 19:05:06 garyb Exp $ */

#ifndef lint
static char vcid[] = "$Id: hststuff3.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* hststuff3.c:  
 * 
 * 
 */
#include "orbfit.h"

/*suffix for orbit-fit files*/
#define SUFFIX ".abg"
#define OBSCODE_HST 2000
#define TIME_STEP   0.2	  /*time step (days) to locate oppangle*/
#define MAX_EPOCHS  100
#define LAT_MAX     15     /*only objects that are this close to ecliptic*/
#define OBSCODE_GRND 807	/*Where the followup observation will be*/
#define ERROR_GRND  0.2	   /*Astrometric error on followup*/
#define PHASE_STEPS 20	   /*number of phases for each object's orbit*/

char   *help[] = {
    "hststuff3:  show position diffusion of a population of objects."
    "            This time put each object at several points in its orbit."
    " usage: hststuff2 oppangle day1 day2 ... dayn",
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
  PBASIS p;
  ORBIT o;
  XVBASIS xv;

  OBSERVATION	futobs;
  FILE *abgfile;
  char	inbuff[256],objname[256],fname[256];

  double **covar,**sigxy;
  double lat,lon;
  int i,nfields;
  double x[3], v[3], range;

  double oppangle, epochs[MAX_EPOCHS], observe_day, period;
  double x0, y0;
  int nepochs, iphase;

  if (argc<3) print_help();
  oppangle = atof(argv[1]);
  for (i=2; i<argc; i++) epochs[i-2] = atof(argv[i]);
  nepochs = argc-2;

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

  printf("#   ID    epoch     x      y     a      e     i    d \n");

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

    for (iphase=0; iphase<PHASE_STEPS; iphase++) {
      /* Advance the object through its orbit by moving back time
       * of perihelion (which is in days) */
      o.T -= period/ (PHASE_STEPS*DAY);
      /* Set up new orbit system based on this */
      elements_to_pbasis(&o, jd0, futobs.obscode, &p);

      
      /* Now look for the day at which opp. angle meets spec */
      observe_day = jd0;

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
    
      x0 = futobs.thetax; y0=futobs.thetay;

    /* print out position at all epochs */
      for (i=0; i<nepochs; i++) {
	futobs.obstime = (observe_day + epochs[i]-jd0)*DAY;
	futobs.xe = -999.;		/* Force evaluation of earth3d */
	predict_posn(&p,NULL,&futobs,NULL);
      
	/* Now transform ecliptic*/
	proj_to_ec(futobs.thetax,futobs.thetay,
		   &lat, &lon,
		   lat0, lon0, NULL);

	/* Get the distance from Earth at observation */
	kbo3d(&p, futobs.obstime, x, v, NULL, NULL, NULL);
	range = (x[0]-futobs.xe)*(x[0]-futobs.xe) +
	  (x[1]-futobs.ye)*(x[1]-futobs.ye) +
	  (x[2]-futobs.ze)*(x[2]-futobs.ze);
	range = sqrt(range);

	printf("%8s %2d %5.1f %7.1f %7.1f %4.1f %5.3f %5.0f %4.1f\n",
	       objname, iphase, epochs[i],
	       (futobs.thetax - x0)/ARCSEC,
	       (futobs.thetay - y0)/ARCSEC,
	       o.a, o.e, o.i, range);
      } /*epoch loop*/
    } /*phase loop */
  } /*object loop*/
  exit(0);
}
