/* 	$Id: mrqcof_orbit.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */
#ifndef lint
static char vcid[] = "$Id: mrqcof_orbit.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */

/* Altered version of Numerical Recipes "mrqcof" to deal with our situation
 * in which we are simultaneously fitting the x and y data.  Also have
 * made use of the existing "OBSERVATION" structure as input data.
 * 4/14/99 gmb
 * 6/13/00 gmb:  add optional extra "observation" that pushes solutions
 * to be bound & nearly circular when close to degenerate.
 * see CVS logs for further alteration comments.
 */
/* energy_wt gives the relative weighting to the binding constraint.  
 * at energy_wt=1., unbound orbit is considered 2-sigma deviation from
 * the circular ideal (as is an aphelic plunging orbit.)
 * energy_wt=0 ignores this, value >1 is stronger incentive to circular.
 */

#define NRANSI
#include "nrutil.h"
#include "orbfit.h"

void mrqcof_orbit(OBSERVATION obsarray[],
		  int ndata, double a[], int ia[],
		  int ma, double **alpha, double beta[], double *chisq,
		  double energy_wt)
{
	int i,j,k,l,m,mfit=0;
	double xmod,ymod,wtx,wty,sig2x,sig2y,dx,dy,*dyda,*dxda;
	double fb1;
	PBASIS	params;
	OBSERVATION *oo;

	dyda=dvector(1,ma);
	dxda=dvector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;

	params.a = a[1];
	params.adot = a[2];
	params.b = a[3];
	params.bdot = a[4];
	params.g = a[5];
	params.gdot = a[6];

	for (i=1;i<=ndata;i++) {
	        oo = &obsarray[i-1];
		kbo2d(&params,oo,&xmod,dxda,&ymod,dyda);
		sig2x=1.0/(oo->dthetax*oo->dthetax);
		sig2y=1.0/(oo->dthetay*oo->dthetay);
		dx=oo->thetax-xmod;
		dy=oo->thetay-ymod;
		for (j=0,l=1;l<=ma;l++) {
		  if (ia[l]) {
		    wtx=dxda[l]*sig2x;
		    wty=dyda[l]*sig2y;
		    for (j++,k=0,m=1;m<=l;m++)
		      if (ia[m]) alpha[j][++k] += wtx*dxda[m]+wty*dyda[m];
		    beta[j] += dx*wtx+dy*wty;
		  }
		}
		*chisq += dx*dx*sig2x+dy*dy*sig2y;
	}
	/*** Add here the energy terms:  the "data" in this context is
	 * the fractional diff fb of tranverse KE from -PE/2. */
	if (energy_wt>0.) {
	  double g3, g4;
	  g3 = pow(params.g, -3.)/GM;
	  g4 = g3/params.g;
	  fb1 = (params.adot*params.adot 
		 + params.bdot*params.bdot
		 + params.gdot*params.gdot) *  g3
	        -1.;
	  sig2x = 4.*energy_wt;	  /*Make fb1=1 be 2-sigma*/
	  dxda[1] = dxda[3] = 0.;
	  dxda[2] = 2*params.adot * g3;
	  dxda[4] = 2*params.bdot * g3;
	  dxda[6] = 2*params.gdot * g3;
	  dxda[5] = -3.*(params.adot*params.adot 
			 + params.bdot*params.bdot
			 + params.gdot*params.gdot) * g4;
	  for (j=0,l=1;l<=ma;l++) {
		  if (ia[l]) {
		    wtx=dxda[l]*sig2x;
		    for (j++,k=0,m=1;m<=l;m++)
		      if (ia[m]) alpha[j][++k] += wtx*dxda[m];
		    beta[j] += -fb1*wtx;
		  }
	  }
	  *chisq += fb1*fb1*sig2x;

	  /* Add 2nd-derivative terms for the energy constraint */
	  sig2x *= fb1;
	  alpha[2][2] += sig2x * 2.* g3;
	  alpha[4][4] += sig2x * 2.* g3;
	  alpha[6][6] += sig2x * 2.* g3;
	  alpha[2][5] -= 6. * params.adot * g4;
	  alpha[5][2] -= 6. * params.adot * g4;
	  alpha[4][5] -= 6. * params.bdot * g4;
	  alpha[5][4] -= 6. * params.bdot * g4;
	  alpha[6][5] -= 6. * params.gdot * g4;
	  alpha[5][6] -= 6. * params.gdot * g4;
	  alpha[5][5] += 12.*(fb1+1.)/(params.g*params.g);
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	free_dvector(dxda,1,ma);
	free_dvector(dyda,1,ma);
}
#undef NRANSI
