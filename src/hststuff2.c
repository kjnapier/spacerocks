/* 	$Id: hststuff2.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: hststuff2.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* planner.c:  give positions and uncertainties for a list of objects
 * at a given site & time. (derived from predict.c)
 * 6/20/00 gmb
 */
#include "orbfit.h"

/*suffix for orbit-fit files*/
#define SUFFIX ".abg"
#define OBSCODE_HST 2000
#define TIME_STEP   0.2	  /*time step (days) to locate oppangle*/
#define MAX_EPOCHS  100
#define LAT_MAX     15     /*only objects that are this close to ecliptic*/

char   *help[] = {
    "hststuff2:  show position diffusion of a population of objects."
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

  double oppangle, epochs[MAX_EPOCHS];
  double nextopp, observe_day, difflast, diffnext, lastopp;
  double x0, y0;
  int nepochs;

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

    /* Now look for the day at which opp. angle meets spec */
    lastopp = -999.;
    observe_day = jd0;

    do {
      if (lastopp<-500.) {
	futobs.obstime = (observe_day-jd0)*DAY;
	futobs.xe = -999.;		/* Force evaluation of earth3d */
	predict_posn(&p,covar,&futobs,sigxy);
	lastopp = opposition_angle(&futobs)/DTOR;
      }
      observe_day += TIME_STEP;
      futobs.obstime = (observe_day-jd0)*DAY;
      futobs.xe = -999.;		/* Force evaluation of earth3d */
      predict_posn(&p,covar,&futobs,sigxy);
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

    futobs.obstime = (observe_day-jd0)*DAY;
    futobs.xe = -999.;		/* Force evaluation of earth3d */
    predict_posn(&p,covar,&futobs,sigxy);

    /* Now transform ecliptic*/
    proj_to_ec(futobs.thetax,futobs.thetay,
	       &lat, &lon,
	       lat0, lon0, NULL);

    /* Only keep objects close to ecliptic */
    if (fabs(lat) > LAT_MAX*DTOR) continue;
    
    x0 = futobs.thetax; y0=futobs.thetay;
    /* get orbital elements for this one */
    pbasis_to_bary(&p, &xv, NULL);
    orbitElements(&xv, &o);

    /* print out position at all epochs */
    for (i=0; i<nepochs; i++) {
      futobs.obstime = (observe_day + epochs[i]-jd0)*DAY;
      futobs.xe = -999.;		/* Force evaluation of earth3d */
      predict_posn(&p,covar,&futobs,sigxy);
      
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

      printf("%8s %5.1f %7.1f %7.1f %4.1f %5.3f %5.0f %4.1f\n",
	     objname, epochs[i],
	     (futobs.thetax - x0)/ARCSEC,
	     (futobs.thetay - y0)/ARCSEC,
	     o.a, o.e, o.i/DTOR, range);
    } /*epoch loop*/
  } /*object loop*/
  exit(0);
}
