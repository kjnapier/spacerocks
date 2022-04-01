/* 	$Id: postdict.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: postdict.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* postdict.c - Read a file containing a/b/g orbit fit, then
 * say compare observations to predicted positions.
 * 6/14/00 gmb
 */
#include "orbfit.h"

char   *help[] = {
    "postdict: Read a file containing a/b/g orbit fit, then compare",
    "       predicted RA & dec to real observations.",
    " usage: postdict [-m mpc_error] [-j JPL_file] [-o observatory_file]",
    "                 [-v] <abgfile> ",
    "  mpc_error  is arcsec uncertainty to apply to MPC-format",
    "             observations.  Default is 0.2",
    "  JPL_file   is binary ephemeris file.  Default is binEphem.423, or",
    "             a file specified by environment variable ORBIT_EPHEMERIS",
    "  observatory_file   is file with observatory site data.  Default is",
    "             observatories.dat, or a file specified by environment",
    "             variable ORBIT_OBSERVATORIES",
    "  abgfile    is name of file with orbit info (from fit_radec)",
    "  stdin      are observations in MPC or my format.",
    "  stdout     residuals plus chi-squared values",
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
  OBSERVATION	futobs, obs;
  struct date_time dt;
  char	inbuff[256],rastring[20],decstring[20];
  double **covar,**sigxy,a,b,PA,**derivs;
  double lat,lon,**covecl;
  double ra,dec, **coveq;
  double dx, dy, chisq, elat, elon;
  double xx,yy,xy,bovasqrd,det;
  int i,nfields;

  int iarg=1;
  if (argc>1 && *argv[1]=='^') print_help();
  if (read_options(&iarg, argc, argv)) print_help();
  if (argc-iarg!=1) print_help();
  

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
  derivs = dmatrix(1,2,1,2);
  covar = dmatrix(1,6,1,6);
  covecl = dmatrix(1,2,1,2);
  coveq = dmatrix(1,2,1,2);

  if (read_abg(argv[iarg],&p,covar)) {
    fprintf(stderr, "Error input alpha/beta/gamma file %s\n",argv[iarg]);
    exit(1);
  }

  while (fgets_nocomment(inbuff, 255, stdin, stdout)!=NULL) {
    if (scan_observation(inbuff, &obs)) exit(1);
    /* Set time to years after jd0, rotate to tangent plane coords */
    obs.obstime = (obs.obstime-jd0)*DAY;
    eq_to_ec(obs.thetax,obs.thetay,&elat,&elon,NULL);
    ec_to_proj(elat,elon,&(obs.thetax),&(obs.thetay),
	       lat0,lon0,NULL);
    /* Calculate the position of Earth at this time to avoid doing
     * it many times later: */
    earth3d(obs.obstime, obs.obscode,
	    &(obs.xe),&(obs.ye),&(obs.ze));

    futobs.obstime = obs.obstime;
    futobs.obscode = obs.obscode;
    futobs.xe = -999.;

    predict_posn(&p,covar,&futobs,sigxy);

    /* Errors: */
    dx = obs.thetax - futobs.thetax;
    dy = obs.thetay - futobs.thetay;
    /* Add observational errors into the uncertainty */
    sigxy[1][1] += pow(obs.dthetax,2.);
    sigxy[2][2] += pow(obs.dthetay,2.);
   
    chisq = dx*dx*sigxy[2][2] -2.*dx*dy*sigxy[1][2] + dy*dy*sigxy[1][1];
    chisq /= sigxy[1][1]*sigxy[2][2] - sigxy[1][2]*sigxy[2][1];

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

    printf("%7.4f %9.2f %9.2f %8.3f %8.3f %8.2f %8.2f %7.2f %.2f\n",
	   futobs.obstime, futobs.thetax/ARCSEC, futobs.thetay/ARCSEC,
	   dx/ARCSEC, dy/ARCSEC,
	   a/ARCSEC,b/ARCSEC,PA, chisq
); 

  }
  exit(0);
}
