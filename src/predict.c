/* 	$Id: predict.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: predict.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* predict.c - Read a file containing a/b/g orbit fit, and spit out
 *  predicted RA & dec plus uncertainties on arbitrary date.
 * 8/12/99 gmb 
 */
#include "orbfit.h"

char   *help[] = {
    "predict: Read a file containing a/b/g orbit fit, and spit out",
    "       predicted RA & dec plus uncertainties on arbitrary date.",
    " usage: predict  [-j JPL_file] [-o observatory_file] [-v] <abgfile>",
    "  JPL_file   is binary ephemeris file.  Default is binEphem.423, or",
    "             a file specified by environment variable ORBIT_EPHEMERIS",
    "  observatory_file   is file with observatory site data.  Default is",
    "             observatories.dat, or a file specified by environment",
    "             variable ORBIT_OBSERVATORIES",
    "  abgfile    is name of file with orbit info (from fit_radec)",
    "  stdin      is list of times to predict, one per line.  Times",
    "             may be given as JD, or YYYY MM DD[.DDDD] [HH MM SS.SS]",
    "             First input line gives observatory code.",
    "  stdout     is predicted RA & Dec plus error ellipse.",
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

  int iarg=1;
  if (argc>1 && *argv[1]=='^') print_help();
  if (read_options(&iarg, argc, argv)) print_help();
  if (iarg>=argc) print_help();

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

  /* get observatory code */
  fprintf (stderr,"Enter observatory code:\n");
  if (fgets_nocomment(inbuff,255,stdin,NULL)==NULL
      || sscanf(inbuff,"%d",&futobs.obscode)!=1) {
    fprintf(stderr,"Error reading observatory code\n");
    exit(1);
  }
  printf("# For observations at site %d\n"
	 "#                 x/RA           y/DEC          "
	 "err_a    err_b  err_pa\n",futobs.obscode);

  fprintf (stderr,"Enter JD's or Y M D ... of observations, -1 to quit:\n");
  while ( fgets_nocomment(inbuff,255,stdin,NULL)!=NULL) {
    nfields=sscanf(inbuff,"%lf %lf %lf %lf %lf %lf",
		   &yr,&mo,&day,&hr,&mn,&ss);
    if (nfields==0 ) {
      fprintf(stderr,"Error on time spec:\n->%s\n",inbuff);
      exit(1);
    } else if (yr<0.) {
      /*done*/
      exit(0);
    } else if (nfields==1 || nfields==2) {
      /* Got a JD. (probably...)*/
      futobs.obstime = (yr-jd0)*DAY;
    } else {
      dt.y = yr;
      dt.mo = mo;
      dt.d = day;
      if (nfields>=4) dt.h = hr;  else dt.h=0.;
      if (nfields>=5) dt.mn = mn; else dt.mn=0.;
      if (nfields>=6) dt.s = ss;  else dt.s=0.;
      futobs.obstime = (date_to_jd(dt)-jd0)*DAY;
    }

    futobs.xe = -999.;		/* Force evaluation of earth3d */

    printf("At time= %s",inbuff);

    predict_posn(&p,covar,&futobs,sigxy);

    printf("# Solar Elongation = %.2f Opp angle = %.2f\n",
	   elongation(&futobs)/DTOR,opposition_angle(&futobs)/DTOR);

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

    printf("Cum. ecliptic posn: %10.4f %10.4f   %8.2f %8.2f %7.2f\n",
	   futobs.thetax/ARCSEC, futobs.thetay/ARCSEC,
	   a/ARCSEC,b/ARCSEC,PA); 

    /* Now transform to RA/DEC, via ecliptic*/
    proj_to_ec(futobs.thetax,futobs.thetay,
	       &lat, &lon,
	       lat0, lon0, derivs);
    /* map the covariance */
    covar_map(sigxy, derivs, covecl, 2, 2);

    /* Now to ICRS: */
    ec_to_eq(lat, lon, &ra, &dec, derivs);
    /* map the covariance */
    covar_map(covecl, derivs, coveq, 2, 2);

    /* Compute a, b, theta of error ellipse for output */
    xx = coveq[1][1]*cos(dec)*cos(dec);
    xy = coveq[1][2]*cos(dec);
    yy = coveq[2][2];
    PA = 0.5 * atan2(2.*xy,(xx-yy)) * 180./PI;	/*go right to degrees*/
    /* Put PA N through E */
    PA = 90.-PA;
    bovasqrd  = (xx+yy-sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) 
      / (xx+yy+sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) ;
    det = xx*yy-xy*xy;
    b = pow(det*bovasqrd,0.25);
    a = pow(det/bovasqrd,0.25);

    ra /= DTOR;
    if (ra<0.) ra+= 360.;
    dec /= DTOR;
    deghms(ra,rastring);
    degdms(dec,decstring);
    printf("ICRS position: %s %s %10.4f %10.4f %7.2f\n",
	   rastring,decstring,a/ARCSEC,b/ARCSEC,PA); 
  }
  exit(0);
}
