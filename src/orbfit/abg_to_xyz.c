/* 	$Id: abg_to_xyz.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: abg_to_xyz.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* abg_to_xyz.c - program to tranform an orbit fit in alpha,beta,
   gamma basis (and its covariance matrix) into barycentric ecliptic
   Cartesian basis of x,y,z,xdot,ydot,zdot.
   8/12/99 gmb
*/

#include "orbfit.h"
char   *abg_to_xyz_help[] = {
    "abg_to_xyz: tranforms an orbit fit in alpha,beta, gamma basis ",
    "         (and its covariance matrix) into barycentric ecliptic",
    "         Cartesian basis of x,y,z,xdot,ydot,zdot",
    " usage:  <abgfile abg_to_xyz >xyzfile",
    0
};

int
abg_to_xyz(int argc, char *argv[])
{
  PBASIS p;
  XVBASIS xv;
  double  **covar_abg, **covar_xyz, **derivs;

  int	i,j;

/*  if (argc!=1) print_help(abg_to_xyz_help); */

  covar_abg = dmatrix(1,6,1,6);
  covar_xyz = dmatrix(1,6,1,6);
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

  /* Print out the results, with comments */
  printf("#     x              y             z           "
	 "vx          vy            vz\n");
  printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n",
	 xv.x, xv.y, xv.z, xv.xdot, xv.ydot, xv.zdot);
  printf("# covariance matrix:\n");
  print_matrix(stdout,covar_xyz,6,6);
  printf("# Position is for JD:\n%.5f\n",jd0);

  free_dmatrix(covar_abg,1,6,1,6);
  free_dmatrix(covar_xyz,1,6,1,6);
  free_dmatrix(derivs,1,6,1,6);
  exit(0);
}
