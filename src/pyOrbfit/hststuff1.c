/* hststuff1.c - Read a file containing a/b/g orbit fit, then
 * spit out a set of mock HST observations at specified opposition
 * angles.
 */
#include "orbfit.h"

#define OBSCODE_HST 2000
/* Start observations today (8/21/01) */
#define JDSTART 2452142.
#define ZA_LIMIT (90.*DTOR)	/*maximum zenith angle for observations*/
#define NORBITS  10	/*number of orbits to observe*/
#define NFIELDS  4      /*number of fields to observe in cycle*/
#define EXP_PER_ORBIT 6	/*exposures per orbit*/
#define ASTROM_ERROR  0.004	/*uncertainty per visit */
#define TIME_STEP (1./24./60.)	/*1-minute time steps to check visibility*/

char   *help[] = {
  "hststuff1: pretend HST observations.",
  "usage:  hststuff1 opp1 opp2 ... oppN",
  "input:  abg file for object",
  "output: positions for HST observations taken at each opp angle",
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
  OBSERVATION	futobs;
  struct date_time dt;
  char	inbuff[256],rastring[20],decstring[20];
  double **covar,**sigxy,a,b,PA,**derivs;
  double lat,lon,**covecl;
  double ra,dec, **coveq;
  double yr,mo,day,hr,mn,ss;
  double xx,yy,xy,bovasqrd,det;
  int i,nfields;

  double observe_day, vis_start, vis_end, jd;
  double lastopp = -999., diffnext, difflast;
  int iorbit, iexpose, x0=-999., y0=-999.;

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

  covar = dmatrix(1,6,1,6);
  sigxy = dmatrix(1,2,1,2);
  if (read_abg(NULL,&p,covar)) {
    fprintf(stderr, "Error in input alpha/beta/gamma file\n");
    exit(1);
  }

  /* get observatory code */
  /** stick with HST: **/
  futobs.obscode = OBSCODE_HST;
  observe_day = JDSTART;

  for (i=1; i<argc; i++) {
    double oppangle, nextopp;
    /* get the next desired opp angle, and advance to the day on which
     * it occurs. */
    if ( sscanf(argv[i],"%lf",&oppangle)!=1) print_help();
    /* Now look for the day at which opp. angle meets spec */
    do {
      if (lastopp<-500.) {
	futobs.obstime = (observe_day-jd0)*DAY;
	futobs.xe = -999.;		/* Force evaluation of earth3d */
	predict_posn(&p,covar,&futobs,sigxy);
	lastopp = opposition_angle(&futobs)/DTOR;
      }
      observe_day += 1.;
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

    /* Advance time until object is NOT visible */
    vis_end = observe_day;
    do {
      vis_end+= TIME_STEP;
      futobs.obstime = (vis_end-jd0)*DAY;
      futobs.xe = -999.;		/* Force evaluation of earth3d */
      predict_posn(&p,covar,&futobs,sigxy);
    } while (zenith_angle(&futobs) <= ZA_LIMIT) ;

    for (iorbit=0; iorbit<NORBITS*NFIELDS; iorbit++) {
      /* Then advance to first time it IS visible. */
      vis_start = vis_end;
      do {
	vis_start+= TIME_STEP;
	futobs.obstime = (vis_start-jd0)*DAY;
	futobs.xe = -999.;		/* Force evaluation of earth3d */
	predict_posn(&p,covar,&futobs,sigxy);
      } while (zenith_angle(&futobs) > ZA_LIMIT) ;

      /* Now find end of visibility period */
      vis_end = vis_start;
      do {
	vis_end+= TIME_STEP;
	futobs.obstime = (vis_end-jd0)*DAY;
	futobs.xe = -999.;		/* Force evaluation of earth3d */
	predict_posn(&p,covar,&futobs,sigxy);
      } while (zenith_angle(&futobs) <= ZA_LIMIT) ;
      vis_end -= TIME_STEP;

      /* Only proceed if we are observing this field on this orbit*/
      if (iorbit%NFIELDS !=0) continue;

      printf("# Visibility on orbit %d is %.2f minutes\n", iorbit,
	     (vis_end - vis_start)*24.*60.);

      /* Find middle of each exposure time */
      jd = vis_start + (vis_end - vis_start)/(2.*EXP_PER_ORBIT);
      for (iexpose=0; iexpose<EXP_PER_ORBIT; iexpose++, 
	     jd+=(vis_end - vis_start)/EXP_PER_ORBIT) {
	futobs.obstime = (jd-jd0)*DAY;
	futobs.xe = -999.;		/* Force evaluation of earth3d */
	predict_posn(&p,covar,&futobs,sigxy);
	
	/* Now transform to RA/DEC, via ecliptic*/
	proj_to_ec(futobs.thetax,futobs.thetay,
		   &lat, &lon,
		   lat0, lon0, NULL);
	/* Now to ICRS: */
	ec_to_eq(lat, lon, &ra, &dec, NULL);
	ra /= DTOR;
	if (ra<0.) ra+= 360.;
	dec /= DTOR;
	deghms(ra,rastring);
	degdms(dec,decstring);
	printf("%14.6f %16s %16s %9.4f %5d\n",
	       jd,rastring,decstring,
	       ASTROM_ERROR*sqrt(NORBITS*EXP_PER_ORBIT),
	       OBSCODE_HST);
	/* Print other useful information*/
	if (x0==-999.) x0=futobs.thetax/ARCSEC;
	if (y0==-999.) y0=futobs.thetay/ARCSEC;
	printf("# posn %10.4f %10.4f Elong %6.2f Opp %7.2f ZA %6.2f\n",
	       futobs.thetax/ARCSEC-x0, futobs.thetay/ARCSEC-y0,
	       elongation(&futobs)/DTOR,
	       opposition_angle(&futobs)/DTOR,
	       zenith_angle(&futobs)/DTOR);
      } /*exposure loop*/
    } /*orbit loop*/
  } /*opposition-angle loop*/
  exit(0);
}
