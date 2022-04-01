/* 	$Id: orbfit2.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: orbfit2.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/***** Routines that execute the fit of orbit to observation data. ****/
/*** assembled 6/14/00 gmb from subroutines scattered about elsewhere.
** will use the mrqmin_orbit and mrqcof_orbit subroutines as well.
**/
#include "orbfit.h"

#define  MAXIT	100	/*maximum mrqmin iterations*/
#define	 CHITOL 0.0001	/*tolerance for chisq minimum*/

/*If this many sigma^2 of gdot is unbound orbit, use energy constraint*/
#define  SIGMA2_USE_ENERGY 9.

/*if fractional error in linear-fit gamma exceeds this,
 *set gdot=0 to avoid degenerate fit: */
#define  GDEGENERATE  0.2  

/*For gdot=0, how many sigma^2 is assigned to unbound gdot?*/
#define  SIGMA2UNBOUND 3.

void 
mrqmin_orbit(OBSERVATION obsarray[], int ndata, double a[], int ia[],
		  int ma, double **covar, double **alpha, double *chisq,
		  double *alamda, double energy_wt);

void
mrqfit(OBSERVATION *obsarray,
       int nobs,
       PBASIS *p,
       int *ia,
       double **covar,
       double *chisq, 
       double energy_wt);

/* Subroutine which goes through all the steps of fitting a set of
 * observations with an orbit.
 * covar should be 6x6 dmatrix upon entry, will hold covariance on exit.
 * Info on the fit is sent to logfile if it is non-NULL.
 */
/* Return value is 5 if energy constraint is used, 6 otherwise */
int
fit_observations(OBSERVATION obsarray[],
		 int nobs,
		 PBASIS *p,
		 double **covar,
		 double *chisq,
		 int *dof,
		 FILE *logfile)
{
  int	*ia,i,j;
  double *a;
  double gbind2;
  double energy_fit;	/*weight given to binding-energy constraint*/
  int   fitparms=6;
  a     = dvector(1,6);
  ia    = ivector(1,6);
  energy_fit = 0.;

  /* Don't even try the no-energy fit unless there are 3 observations*/
  if (nobs<2) {
    fprintf(stderr,"ERROR: not enough observations nobs=%d\n",nobs);
  } else if (nobs==2) {
    fitparms=4;
  } else {
    prelim_fit(obsarray,nobs,p,covar);
    /* If the preliminary linear fit has left the distance indeterminate
     * with the simple 5-dimensional solution, then this data is probably
     * doubly indeterminate and we must jump to the gdot=0 + energy 
     * constraint.
     */
    if (covar[5][5]<0. || covar[5][5]/(p->g*p->g)>GDEGENERATE*GDEGENERATE) 
      fitparms = 4;
  }

  if (fitparms==6) {
    if (logfile!=NULL) {
      fprintf(logfile,"#  Preliminary a, adot, b, bdot, g, gdot:\n");
      fprintf(logfile,"#  %lf %lf %lf %lf %lf %lf\n",p->a,p->adot,p->b,
	    p->bdot, p->g, p->gdot);  
    }
    /* Prepare for the first Marquandt minimization */
    for (i=1; i<=6; i++) ia[i]=1;
    
    mrqfit(obsarray, nobs, p, ia, covar, chisq, energy_fit);
    
    /* See if the covar matrix allows gdot to be large enough
     * to allow unbound orbits.  If so, refit with energy constraint.
     */
    {
      double t, var_gbind2;
      /* this is the square of the binding limit on gdot:*/
      t = xBary*xBary + yBary*yBary + zBary*zBary;
      t = 1 + t*p->g*p->g - 2*p->g*zBary; /*t=g^2 * barycentric dist^2 */
      gbind2 = 2*GM*pow(p->g,3.)/sqrt(t) - p->adot*p->adot - p->bdot*p->bdot;
      /* Allow a 2-sigma excursion: */
      var_gbind2 = 36*GM*GM*pow(p->g,4.)/t*covar[5][5] +
	4*p->adot*p->adot*covar[2][2] + 
	4*p->bdot*p->bdot*covar[4][4];
      gbind2 -= 2*sqrt(var_gbind2);	/*take a lower limit here*/
    }


    if (SIGMA2_USE_ENERGY*covar[6][6] >= gbind2) {
      /* If unconstrained gdot variance allows unbound,
       * fit again with energy constrained.
       */
      fitparms=5;
    }
  }

  /* Execute the energy-constrained fit if nobs==2 or the unconstrained
   * fit demonstrated degeneracy in gdot: */
  if (fitparms<6) {

    if (logfile!=NULL) 
      fprintf(logfile, "#WARNING: Fitting with energy constraint\n");

    /* Prepare for the first Marquandt minimization */
    for (i=1; i<=6; i++) ia[i]=1;

    if (fitparms==4) {
      /* avoid the (degenerate) preliminary fit, just choose values*/
      p->g = 0.03;
      p->a = p->b = p->bdot = p->gdot = 0.;
      p->adot = 0.03;
      /* And don't fit gdot at all: */
      ia[6] = 0;
      if (logfile!=NULL)
	fprintf(logfile,"#WARNING:  and gdot fixed =0\n");
      
    } else {
      prelim_fit(obsarray,nobs,p,covar);
    }

    /* Fit first with strong weight to circular: */
    energy_fit = 10.;
    mrqfit(obsarray, nobs, p, ia, covar, chisq, energy_fit);
    /* Now relax to get a good error estimate*/
    energy_fit = 1.;
    mrqfit(obsarray, nobs, p, ia, covar, chisq, energy_fit);
  }

  free_dvector(a,1,6);
  free_ivector(ia,1,6);

  if (fitparms==4) {
    /* put nominal uncertainties on gdot */
    covar[6][6] = GM*pow(p->g,3.)/SIGMA2UNBOUND;
     for (i=1; i<=5; i++) 
       covar[i][6]=covar[6][i]=0.;
  }

  *dof = 2*nobs - fitparms;
  return fitparms;
}


/* Subroutine which executes the mrqmin optimization */
void
mrqfit(OBSERVATION *obsarray,
       int nobs,
       PBASIS *p,
       int *ia,
       double **covar,
       double *chisq,
       double energy_wt)
{
   double *a, alambda,oldchi,**alpha;
   int  ma=6,niter;

   a     = dvector(1,6);
   alpha = dmatrix(1,6,1,6);

   a[1]=p->a;
   a[2]=p->adot;
   a[3]=p->b;
   a[4]=p->bdot;
   a[5]=p->g;
   a[6]=p->gdot;

   alambda=-1.;
   *chisq = 1.e14;
   niter = 0;
   do {
     /*do an iteration of the nr marquandt search*/
     oldchi = *chisq;
     niter++;
     mrqmin_orbit(obsarray, nobs, a, ia,
		  ma, covar, alpha, chisq, &alambda, energy_wt);
#ifdef DEBUG
     fprintf(stderr,
	     "# Iteration %d: chisq %f, alamda now %f\n",
	     niter,*chisq,alambda);  
     fprintf(stderr,"#  Params %g %g %g %g %g %g\n",
		 a[1],a[2],a[3],a[4],a[5],a[6]);
#endif
     if (alambda>1e8) {
      /*
       fprintf(stderr,
	       "MRQMIN being stopped after %d iterations, alambda=%f\n",
	       niter,alambda);
      */
       break;
     }
   } while (niter<MAXIT && (*chisq<oldchi-CHITOL || *chisq>=oldchi));
   /*first get the uncertainties back from mrqmin*/
   alambda=0.;
   mrqmin_orbit(obsarray, nobs, a, ia,
		ma, covar, alpha, chisq, &alambda, energy_wt);

   p->a = a[1];
   p->adot = a[2];
   p->b = a[3];
   p->bdot = a[4];
   p->g = a[5];
   p->gdot = a[6];

   return;
}


