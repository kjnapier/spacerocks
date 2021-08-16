/* Version of Numerical Recipes mrqmin.c altered to deal specifically
 * with the orbit fitting.  The only real change here is in the call
 * to mrqcof, and the elimination of a number of arrays that are passed
 * through to mrqcof in favor of our OBSERVATIONS structures.
 *  Also all floats are now doubles.

 * Alteration 6/13/00 gmb:  add a flag on input which will tell it to
 * add an additional constraint (disguised as an observation) that
 * the orbit be near circular or at least bound.
 */
#define NRANSI
#include "nrutil.h"
#include "orbfit.h"

void mrqmin_orbit(OBSERVATION obsarray[], int ndata, double a[], int ia[],
		  int ma, double **covar, double **alpha, double *chisq,
		  double *alamda, double energy_wt)
{
	void covsrt(double **covar, int ma, int ia[], int mfit);
	void gaussj(double **a, int n, double **b, int m);
	void mrqcof_orbit(OBSERVATION obsarray[],
			  int ndata, double a[], int ia[],
			  int ma, double **alpha, double beta[], double *chisq,
			  double energy_wt);
	int j,k,l,m;
	static int mfit;
	static double ochisq,*atry,*beta,*da,**oneda;

	if (*alamda < 0.0) {
		atry=dvector(1,ma);
		beta=dvector(1,ma);
		da=dvector(1,ma);
		for (mfit=0,j=1;j<=ma;j++)
			if (ia[j]) mfit++;
		oneda=dmatrix(1,mfit,1,1);
		*alamda=0.001;
		mrqcof_orbit(obsarray,ndata,a,ia,ma,alpha,beta,chisq,
			     energy_wt);
		ochisq=(*chisq);
		for (j=1;j<=ma;j++) atry[j]=a[j];
	}
	for (j=0,l=1;l<=ma;l++) {
		if (ia[l]) {
			for (j++,k=0,m=1;m<=ma;m++) {
				if (ia[m]) {
					k++;
					covar[j][k]=alpha[j][k];
				}
			}
			covar[j][j]=alpha[j][j]*(1.0+(*alamda));
			covsrt(alpha,ma,ia,mfit); /**/
			oneda[j][1]=beta[j];
		}
	}
	gaussj(covar,mfit,oneda,1);
	for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0) {
		covsrt(covar,ma,ia,mfit);
		free_dmatrix(oneda,1,mfit,1,1);
		free_dvector(da,1,ma);
		free_dvector(beta,1,ma);
		free_dvector(atry,1,ma);
		return;
	}
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) atry[l]=a[l]+da[++j];
		mrqcof_orbit(obsarray,ndata,atry,ia,ma,covar,da,chisq,
			     energy_wt);
	if (*chisq < ochisq) {
		*alamda *= 0.1;
		ochisq=(*chisq);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				for (j++,k=0,m=1;m<=ma;m++) {
					if (ia[m]) {
						k++;
						alpha[j][k]=covar[j][k];
					}
				}
				beta[j]=da[j];
				a[l]=atry[l];
			}
		}
	} else {
		*alamda *= 10.0;
		*chisq=ochisq;
	}
}
#undef NRANSI
