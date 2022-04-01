/* 	$Id: abg_to_aei.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: abg_to_aei.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* abg_to_aei.c - program to tranform an orbit fit in alpha,beta,
   gamma basis (and its covariance matrix) into barycentric ecliptic
   Cartesian basis of x,y,z,xdot,ydot,zdot.
   8/12/99 gmb
   **** and thence into true orbital parameters.  10/20/99 gmb
*/

#include "orbfit.h"
char   *abg_to_aei_help[] = {
    "abg_to_aei: tranforms an orbit fit in alpha,beta, gamma basis ",
    "         (and its covariance matrix) into usual orbital parameters",
    " usage:  <abgfile abg_to_xyz >aeifile",
    0
};

int
abg_to_aei(int argc, char *argv[])
{
  PBASIS p;
  XVBASIS xv;
  ORBIT orbit;
  double  **covar_abg, **covar_xyz, **derivs, **covar_aei;

  int	i,j;

/*  if (argc!=1) print_help(abg_to_aei_help); */

  covar_abg = dmatrix(1,6,1,6);
  covar_xyz = dmatrix(1,6,1,6);
  covar_aei = dmatrix(1,6,1,6);
  derivs = dmatrix(1,6,1,6);

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

  if (read_abg(NULL,&p,covar_abg)) {
    fprintf(stderr, "Error in input alpha/beta/gamma file\n");
    exit(1);
  }

  /* Transform the orbit basis and get the deriv. matrix */
  pbasis_to_bary(&p, &xv, derivs);

  /* Map the covariance matrix to new basis */
  covar_map(covar_abg, derivs, covar_xyz,6,6);

  /* Get partial derivative matrix from xyz to aei */
  aei_derivs(&xv, derivs);
  /* Map the covariance matrix to new basis */
  covar_map(covar_xyz, derivs, covar_aei,6,6);

  /* Transform xyz basis to orbital parameters */
  orbitElements(&xv, &orbit);


  /* Print out the results, with comments */
  printf("# Barycentric osculating elements in ICRS at epoch %.1f:\n",jd0);
  printf("#    a            e       i      Node   Arg of Peri   Time of Peri\n");
  printf("%12.6f  %9.6f  %8.3f %8.3f  %8.3f %11.3f\n",
	 orbit.a, orbit.e, orbit.i, orbit.lan, orbit.aop, orbit.T);
  printf("+-%10.6f  %9.6f  %8.3f %8.3f  %8.3f %11.3f\n",
	 sqrt(covar_aei[1][1]),
	 sqrt(covar_aei[2][2]),
	 sqrt(covar_aei[3][3])/DTOR,
	 sqrt(covar_aei[4][4])/DTOR,
	 sqrt(covar_aei[5][5])/DTOR,
	 sqrt(covar_aei[6][6])/DAY);
  printf("# covariance matrix:\n");
  print_matrix(stdout,covar_aei,6,6);

  free_dmatrix(covar_abg,1,6,1,6);
  free_dmatrix(covar_xyz,1,6,1,6);
  free_dmatrix(covar_aei,1,6,1,6);
  free_dmatrix(derivs,1,6,1,6);
  exit(0);
}
