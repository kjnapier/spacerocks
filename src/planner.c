/* 	$Id: planner.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: planner.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* planner.c:  give positions and uncertainties for a list of objects
 * at a given site & time. (derived from predict.c)
 * 6/20/00 gmb
 */
#include <string.h>
#include "orbfit.h"

/*suffix for orbit-fit files*/
#define SUFFIX ".abg"

char   *help[] = {
    "planner:  give positions and uncertainties for a list of objects",
    "       on a given date & site",
    " usage: planner [-j JPL_file] [-o observatory_file] [-v] <siteid> <jd>",
    "              OR ",
    "        planner [-j JPL_file] [-o observatory_file] <siteid> <year> <month> <date>",
    "  JPL_file   is binary ephemeris file.  Default is binEphem.423, or",
    "             a file specified by environment variable ORBIT_EPHEMERIS",
    "  observatory_file   is file with observatory site data.  Default is",
    "             observatories.dat, or a file specified by environment",
    "             variable ORBIT_OBSERVATORIES",
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
  OBSERVATION	futobs;
  FILE *abgfile;
  struct date_time dt;
  char	inbuff[256],objname[256],fname[256];
  char	rastring[32],decstring[32];

  double **covar,**sigxy,a,b,PA,**derivs;
  double lat,lon,**covecl;
  double ra,dec, **coveq;
  double jd;
  double xx,yy,xy,bovasqrd,det;
  int i,nfields, siteid;
  double x[3], v[3], range;
  int iarg=1;
  if (argc>1 && *argv[1]=='^') print_help();
  if (read_options(&iarg, argc, argv)) print_help();
  if (iarg>=argc) print_help();

  siteid = atoi(argv[iarg]);
  if (siteid<=0) {
    fprintf(stderr,"Invalid observing site ID %s\n",argv[iarg]);
    print_help();
  }
  futobs.obscode = siteid;
  iarg++;

  if (argc-iarg==1) {
    jd = atof(argv[iarg]);
  } else if (argc-iarg==3) {
    dt.y = atoi(argv[++iarg]);
    dt.mo = atoi(argv[++iarg]);
    dt.d = atof(argv[iarg]);
    dt.h = dt.mn = dt.s = 0.;
    if (dt.y <= 0 || dt.mo <=0 || dt.mo > 31) {
      fprintf(stderr,"Invalid observing date ID %s\n",argv[1]);
      print_help();
    }
    jd = date_to_jd(dt);
  } else {
    fprintf(stderr,"Incorrect number of arguments\n");
    print_help();
  }

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

  printf("# For observations at site %d at JD %.3f:\n",
	 futobs.obscode, jd);
  printf("#   ID    elong opp range    ra    (J2000)"
	 "   dec          a        b     PA\n");

  sigxy = dmatrix(1,2,1,2);
  derivs = dmatrix(1,2,1,2);
  covar = dmatrix(1,6,1,6);
  covecl = dmatrix(1,2,1,2);
  coveq = dmatrix(1,2,1,2);

  /* Loop through input containing names of objects */
  while ( fgets_nocomment(inbuff,255,stdin,NULL)!=NULL) {
    sscanf(inbuff," %s ",objname);
    strcpy(fname,objname);
    strcat(fname,SUFFIX);

    if (read_abg(fname,&p,covar)) {
      fprintf(stderr, "Error reading alpha/beta/gamma file %s\n",fname);
      exit(1);
    }

    futobs.obstime = (jd-jd0)*DAY;
    futobs.xe = -999.;		/* Force evaluation of earth3d */

    predict_posn(&p,covar,&futobs,sigxy);

    /* Get the distance from Earth at observation */
    kbo3d(&p, futobs.obstime, x, v, NULL, NULL, NULL);
    range = (x[0]-futobs.xe)*(x[0]-futobs.xe) +
      (x[1]-futobs.ye)*(x[1]-futobs.ye) +
      (x[2]-futobs.ze)*(x[2]-futobs.ze);
    range = sqrt(range);


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

    printf("%8s %4.0f %4.0f %4.1f %s %s %8.2f %8.2f %7.2f\n",
	   objname,
	   elongation(&futobs)/DTOR,
	   opposition_angle(&futobs)/DTOR,
	   range,
	   rastring,decstring,a/ARCSEC,b/ARCSEC,PA); 
  }
  exit(0);
}
