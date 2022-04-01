/* 	$Id: AllMPC.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: AllMPC.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* Fit orbit to RA/DEC file.  First
 * preliminary fit, THEN do a full Marquandt-Levenberg fit.
 * Use sigma matrix of full fit to predict positions and
 * errors using the full nonlinear formulae.  
 * usage:  orbtest3 <datafile> <# params>
 *  datafile contains list of observations (e.g. output of orbtest1)
 *  # params is either 5 (don't fit to gdot) or 6 (fit all params)
 * 8/9/99 gmb
 ******
 * This one runs through a file with gobs of objects concatenated in MPC
 * format and fits each one.
 */
#include "orbfit.h"
#include <time.h>
#include <string.h>
/*default V-R color:*/
#define VRCOLOR 0.5

char   *help[] = {
    " AllMPC:  Fit KBO orbits to a long list of MPC observations",
    " usage:  AllMPC [-m mpc_error] [-j JPL_file] [-o observatory_file] [-v]",
    "  mpc_error  is arcsec uncertainty to apply to MPC-format",
    "             observations.  Default is 0.2",
    "  JPL_file   is binary ephemeris file.  Default is binEphem.405, or",
    "             a file specified by environment variable ORBIT_EPHEMERIS",
    "  observatory_file   is file with observatory site data.  Default is",
    "             observatories.dat, or a file specified by environment",
    "             variable ORBIT_OBSERVATORIES",
    "  stdin      is long list of MPC-format observations",
    "  stdout     gives name, distance, & orb elements for each object",
    " An output file <objname>.abg is produced for each object in the",
    " file, giving its orbital parameters & covar matrix in abg form.",
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


int
fit_observations(OBSERVATION obsarray[],
		 int nobs,
		 PBASIS *p,
		 double **covar,
		 double *chisq,
		 int *dof,
		 FILE *logfile);

OBSERVATION obsarray[MAXOBS];
int	nobs;

int
main(int argc, char *argv[])
{
  PBASIS p;
  ORBIT  orbit;
  XVBASIS xv;
  OBSERVATION *obs;

  double **covar, **covar_xyz, **covar_aei, **derivs;
  double chisq, dist, ddist, elat, elon;
  int i, dof, atend, nobj;
  char objname[10]="noname", newname[10],fname[14],inbuff[256];
  FILE *abgfile=NULL;
  double magsum, mag, firstobs, lastobs;
  int nmag, fp;
  char band;

  covar = dmatrix(1,6,1,6);
  covar_xyz = dmatrix(1,6,1,6);
  covar_aei = dmatrix(1,6,1,6);
  derivs = dmatrix(1,6,1,6);

  {
    int iarg=1;
    if (argc>1 && *argv[1]=='^') print_help();
    if (read_options(&iarg, argc, argv)) print_help();
  }

  nobs = 0;
  nobj = 0;

  /* echo the command line to each output file*/
  fprintf(stdout,"#");
  for (i=0; i<argc; i++) fprintf(stdout," %s",argv[i]);
  {
    time_t timettt;
    time(&timettt);
    /* note that ctime returns string with newline at end */
    fprintf(stdout,"\n#---%s",ctime(&timettt));
  }
  fprintf(stdout,"# ID chisq nobs   arc Rmag   bary. dist.       a  "
	  "      e      Sig_aa   Sig_ee   Sig_ae\n");

  /* Start reading the input file */
  while (1) {
    atend = (fgets_nocomment(inbuff,255,stdin,NULL)==NULL);
    if (!atend) strncpy(newname,inbuff+6,6);
    newname[6] = 0;
    if (nobj>0 && (atend || strcmp(newname, objname))) {
      /* Come to the end of data for the previous object.
       * Execute the fitting process and report results.
       */
      
      strcpy(fname, objname);
      strcat(fname, ".abg");
      if ( (abgfile=fopen(fname,"w"))==NULL) {
	fprintf(stderr,"Can't open output file %s\n",fname);
	exit(1);
      }

      /* echo the command line to each output file*/
      fprintf(abgfile,"#");
      for (i=0; i<argc; i++) fprintf(abgfile," %s",argv[i]);
      fprintf(abgfile,"\n# For object code %s\n",objname);
      {
	time_t timettt;
	time(&timettt);
	/* note that ctime returns string with newline at end */
	fprintf(abgfile,"#---%s",ctime(&timettt));
      }

      fprintf(abgfile,"# Fitting %d observations\n",nobs);

      /* Call subroutine to do the actual fitting: */
      fp = fit_observations(obsarray, nobs, &p, covar, &chisq, &dof,abgfile);

      fprintf(abgfile,"# Chi-squared of fit: %.2f DOF: %d\n",chisq,dof);
      fprintf(abgfile,"# Exact a, adot, b, bdot, g, gdot:\n");
      fprintf(abgfile,"%11.8f %11.8f %11.8f %11.8f %11.8f %11.8f\n",
	      p.a,p.adot,p.b,p.bdot, p.g, p.gdot);
      pbasis_to_bary(&p, &xv, derivs);

      orbitElements(&xv, &orbit);
      fprintf(abgfile,"# a=%lf AU,e=%lf,i=%lf deg\n",orbit.a, orbit.e, orbit.i);

      dist = sqrt(xBary*xBary + yBary*yBary + pow(zBary-1/p.g,2.));
      ddist = dist*dist*sqrt(covar[5][5]);
      fprintf(abgfile,"# Barycentric distance %.3f+-%.3f\n",dist,ddist);

      /* Print the covariance matrix to stdout */
      fprintf(abgfile,"# Covariance matrix: \n");
      print_matrix(abgfile,covar,6,6);

      /* Print out information on the coordinate system */
      fprintf(abgfile,
	      "#     lat0       lon0       xBary     yBary      zBary   JD0\n");
      fprintf(abgfile,
	      "%11.6f  %11.6f  %9.6f %9.6f %9.6f  %.5f\n",
	      lat0/DTOR,lon0/DTOR,xBary,yBary,zBary,jd0);
      fclose(abgfile);

      /* get the element uncertainties*/
      covar_map(covar, derivs, covar_xyz,6,6);
      aei_derivs(&xv, derivs);
      covar_map(covar_xyz, derivs, covar_aei,6,6);
      /*dump some info about this object to stdout:*/
      if (nmag>0) 
	magsum /= nmag;
      else
	magsum = 0.;
      printf("%s %3d %5.1f %6.3f %5.2f %5.2f +- %5.2f %10.6f %8.6f %.4g %.4g %.4g %1d\n",
	     objname, nobs, chisq, 
	     lastobs - firstobs ,magsum,
	     dist, ddist, orbit.a, orbit.e, covar_aei[1][1], covar_aei[2][2],
	     covar_aei[1][2], fp);

      /* If that was the last record of the file, we are done*/
      if (atend) exit(0);

      /*Reset the system for a new object */
      nobs = 0;
    }

    /* Add this observation to the list for this object */
    obs = &(obsarray[nobs]);
    if (scan_observation(inbuff, obs)) exit(1);
    eq_to_ec(obs->thetax,obs->thetay,&elat,&elon,NULL);

    if (nobs==0) {
      double xec, yec, zec;
      /* save the name: */
      nobj++;
      strcpy(objname, newname);
      magsum = 0.;
      nmag = 0;

      /* Use first observation to set the reference frame */
      jd0 = obs->obstime;
      lat0 = elat;
      lon0 = elon;
      
      /* Find location of SSBARY wrt observatory at zero time */
      earth_ssbary(jd0, obs->obscode, &xBary, &yBary, &zBary);

      /* Negate the vector to make it earth->SSBARY*/
      /* And rotate equatorial into the tangent-point coords */
      xBary *= -1.;  yBary *= -1.;  zBary *= -1.;
      xyz_eq_to_ec(xBary, yBary, zBary, &xec, &yec, &zec,NULL);
      xyz_ec_to_proj(xec, yec, zec, &xBary, &yBary, &zBary, lat0, lon0, NULL);

      /* reset first/last obs times */
      firstobs=1e9;
      lastobs=-1e9;
    }

    /* Set time to years after jd0, rotate to tangent plane coords */
    obs->obstime = (obs->obstime-jd0)*DAY;
    ec_to_proj(elat,elon,&(obs->thetax),&(obs->thetay),
	       lat0,lon0,NULL);
    /* Calculate the position of Earth at this time to avoid doing
     * it many times later: */
    earth3d(obs->obstime, obs->obscode,
	    &(obs->xe),&(obs->ye),&(obs->ze));
    

    /*Keep track of the mean mag if there are any*/
    if (sscanf(inbuff+65,"%lf %c",&mag,&band)==2) {
      if (band=='R') {
	magsum += mag;
	nmag++;
      } else if (band=='V') {
	magsum += mag-VRCOLOR;
	nmag++;
      }
      else {
	fprintf(stderr,"# photometry in unknown color: %c\n",band);
      }
    }

    /* Keep track of arc length */
    if (obs->obstime < firstobs) firstobs = obs->obstime;
    if (obs->obstime > lastobs) lastobs = obs->obstime;

    nobs++;

  }

  free_dmatrix(covar,1,6,1,6);

  exit(0);
} 

