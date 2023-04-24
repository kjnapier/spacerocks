/* 	$Id: orbfit1.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: orbfit1.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* orbfit1.c - program to fit orbital parameters to observed KBO data.
This is basically a template program which makes some simplifying
assumptions, with the goal of setting up for a more exact fit and for
obtaining estimates of the sizes of error ellipses that we will
encounter in the future.
 
4/9/99 gmb
*/


#include <string.h> 
#include <ctype.h>
#include "orbfit.h"
#include "ephem_types.h"

int invert_matrix(double **in, double **out,int dim);
void lubksb(double **a, int n, int *indx, double b[]);

/* Create the variables which define the coordinate system */
double	lat0, lon0;	/* ecliptic lat & lon of tangent point */
double	xBary, yBary, zBary;	/*Posn of barycenter in our system*/
double	jd0;		/* Zeropoint of time scale */

double mpc_dtheta=DEFAULT_DTHETA;	/*default astrometric error*/

void
set_mpc_dtheta(double d) {
  mpc_dtheta = d;
}

/* Project KBO position onto tangent plane at a given time */
void
kbo2d(PBASIS *pin, 
      OBSERVATION *obs,
      double *x, double dx[],
      double *y, double dy[])
{
  double        xk[3],vk[3],dxk[7],dyk[7],dzk[7];
  double        xe,ye,ze;
  double        invz;
  int           i;
  double distance; 

  /* Get the Earth position if not already calculated*/
  if (obs->xe < -9.) {
    earth3d(obs->obstime, obs->obscode, &xe,&ye,&ze);
    obs->xe = xe;
    obs->ye = ye;
    obs->ze = ze;
  } else {
    xe = obs->xe;
    ye = obs->ye;
    ze = obs->ze;
  }

  /* Get preliminary KBO position */
  kbo3d(pin,obs->obstime,
	xk,vk,dxk,dyk,dzk);

  /* At this point one should account for the light-travel time
     delay between observed time and kbo orbit position time.
     Calculate distance & time delay using naive position.
  */ 
  distance=sqrt( (xk[0]-xe)*(xk[0]-xe) + (xk[1]-ye)*(xk[1]-ye)
		 + (xk[2]-ze)*(xk[2]-ze) );

  kbo3d(pin,obs->obstime-distance/SPEED_OF_LIGHT, 
	xk,vk,dxk,dyk,dzk);
 
  invz = 1./(xk[2]-ze);
  *x = (xk[0] - xe)*invz;
  *y = (xk[1] - ye)*invz;
  /*fprintf(stderr," kbo xyz: %.3f %.3f %.3f earth %.3f %.3f %.3f\n",
    xk[0],xk[1],xk[2],xe,ye,ze);*/
  for (i=1; i<=6; i++) {
    if (dx!=NULL) dx[i] = dxk[i]*invz - (xk[0]-xe)*dzk[i]*invz*invz;
    if (dy!=NULL) dy[i] = dyk[i]*invz - (xk[1]-ye)*dzk[i]*invz*invz;
  }
  /* incorporate derivative w.r.t simplified time delay */
  if (dx!=NULL) dx[5] += vk[0]*invz*invz/SPEED_OF_LIGHT;
  if (dx!=NULL) dy[5] += vk[1]*invz*invz/SPEED_OF_LIGHT;
 
  return;
}
 
/* Linearized version of the 2d projection of orbital position.  Only
   leading-order terms for each parameter's derivative are given.
   The derivative with respect to gdot is inserted here - note that
   it is much smaller than the others, and is the only part of the
   derivatives that has any dependence upon the incoming PBASIS.
*/
void
kbo2d_linear(PBASIS *pin,
	     OBSERVATION *obs,
	     double *x, double dx[],
             double *y, double dy[])
{
  double        xe,ye,ze,t_emit;
  int           i;
 
  /* Get the Earth position if not already calculated*/
  if (obs->xe < -9.) {
    earth3d(obs->obstime, obs->obscode, &xe,&ye,&ze);
    obs->xe = xe;
    obs->ye = ye;
    obs->ze = ze;
  } else {
    xe = obs->xe;
    ye = obs->ye;
    ze = obs->ze;
  }
 
  /* Account for light-travel time differentially to leading order
   * by retarding KBO motion by z component of Earth's position. 
   * Note ignoring acceleration here.
   */
  t_emit = obs->obstime - ze/SPEED_OF_LIGHT;
  *x = pin->a + pin->adot*t_emit - pin->g * xe
    - pin->gdot * (pin->adot*t_emit*t_emit - pin->g*xe*t_emit);
  *y = pin->b + pin->bdot*t_emit - pin->g * ye
    - pin->gdot * (pin->bdot*t_emit*t_emit - pin->g*ye*t_emit);  
 
  if (dx!= NULL && dy!=NULL) {
    for (i=1; i<=6; i++) dx[i] = dy[i]=0.;
 
    dx[1] = dy[3] = 1.;
    dx[2] = dy[4] = t_emit;
    dx[5] = -xe;
    dy[5] = -ye;
    dx[6] = -(pin->adot*t_emit*t_emit - pin->g*xe*t_emit);
    dy[6] = -(pin->bdot*t_emit*t_emit - pin->g*ye*t_emit);
  }
 
  return;
}

/* Read a line from specified input stream, skipping blank or
 * commented (#) lines.  These are echoed to fpout if it's non-null.
 * Returns NULL at EOF.
 */
char *
fgets_nocomment(char *inbuff, int length, 
		FILE *fpin, FILE *fpout)
{
  int needmore=0;
  char test[10];
  while (1) {
    if (fgets(inbuff,length,fpin)==NULL) return(NULL);

    /* Skip blank lines or comments */
    if (sscanf(inbuff,"%8s",test)>0
	&& test[0]!='#') return(inbuff);

    /* Echo comments to output if desired */
    if (fpout!=NULL) 
      fputs(inbuff,fpout);
  }
}

/* read a string as a line of an observation.  Returns jd in the time
 * part of OBSERVATION structure, ra & dec in the x & y parts.  Format
 * error returns non-zero.
 */
int
scan_observation(char *inbuff,
		 OBSERVATION *obs)
{
  char rastring[80],decstring[80],*endptr;
  double jd;
  extern double dmsdeg(char *string);
  extern double hmsdeg(char *string);
  /* get date to see which format this is */
  /* For an MPC format the first field will have non-numeric characters*/
  jd = strtod(inbuff, &endptr);
  if (jd==0. || *endptr==0 || !isspace(*endptr)) {
    /* See if this is perhaps in MPC format */
    struct date_time dd;
    if (sscanf(inbuff+15,"%d %d %f",&(dd.y),&(dd.mo),&(dd.d))!=3) {
      fprintf(stderr,"Format error in observation file:\n ->%s\n",inbuff);
      return(1);
    }
    dd.h = dd.mn = dd.s = 0.;
    jd = date_to_jd(dd);
    strncpy(rastring,inbuff+32,11);
    strncpy(decstring,inbuff+44,11);
    sscanf(inbuff+77,"%3d",&(obs->obscode));
    obs->dthetay = mpc_dtheta;
    
  } else if (jd<10000.) {
    /* See if perhaps this was y/m/d instead of JD: */
    struct date_time dd;
    if (sscanf(inbuff,"%d %d %f %s %s %lf %d",
	       &(dd.y),&(dd.mo),&(dd.d),rastring,decstring,
	       &(obs->dthetay),
	       &(obs->obscode)) !=7 ) {
      fprintf(stderr,"Format error in observation file:\n ->%s\n",inbuff);
      return(1);
    }
    dd.h = dd.mn = dd.s = 0.;
    jd = date_to_jd(dd);
  } else if (sscanf(inbuff,"%lf %s %s %lf %d",
		    &jd,rastring,decstring,
		    &(obs->dthetay),
		    &(obs->obscode)) !=5 ) {
    fprintf(stderr,"Format error in observation file:\n ->%s\n",inbuff);
    return(1);
  }
  
  /* Convert strings to ra & dec */
  obs->obstime = jd;
  obs->thetax = DTOR * hmsdeg(rastring);
  obs->thetay = DTOR * dmsdeg(decstring);
  obs->dthetax = obs->dthetay = obs->dthetay*ARCSEC;
  
  return(0);
}


/* read list of observations from a file */
/* This time we assume that it's in RA & DEC, and we automatically
 * set up the tangent point and jd0 to correspond to 1st observation.
 */
/* Blank & commented lines are skipped */
/* Input file is stdin if fname==NULL */
int
read_radec(OBSERVATION obsarray[], char *fname, int *nobs)
{
  FILE *fptr;
  OBSERVATION  *obs;
  int  i;
  char	inbuff[256];
  double jd,ra,dec,elat,elon;

  if (fname==NULL)
    fptr = stdin;
  else if ( (fptr=fopen(fname,"r"))==NULL) {
    fprintf(stderr,"Error opening observations file %s\n",fname);
    exit(1);
  }

  *nobs=0;
  while ( fgets_nocomment(inbuff,255,fptr,NULL)!=NULL) {
    if ( scan_observation(inbuff, &(obsarray[*nobs]))) {
      fprintf(stderr,"Quitting on format error\n");
      exit(1);
    }

    obs = &(obsarray[*nobs]);
    (*nobs)++;

    eq_to_ec(obs->thetax,obs->thetay,&elat,&elon,NULL);

    if (*nobs==1) {
      double xec, yec, zec;
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
    }

    /* Set time to years after jd0, rotate to tangent plane coords */
    obs->obstime = (obs->obstime-jd0)*DAY;
    ec_to_proj(elat,elon,&(obs->thetax),&(obs->thetay),
	       lat0,lon0,NULL);
    /* Calculate the position of Earth at this time to avoid doing
     * it many times later: */
    earth3d(obs->obstime, obs->obscode,
	    &(obs->xe),&(obs->ye),&(obs->ze));
		
  }
  if (fname!=NULL) fclose(fptr);
  return(0);
}

/* Write an observation (in arcseconds) to already-opened file
 * descriptor.
 */
void
write_obs(FILE *fptr,
	  OBSERVATION *obs)
{
  fprintf(fptr,"%10.7f %9.2f %5.2f %9.2f %5.2f %3d\n",
	  obs->obstime, obs->thetax/ARCSEC, obs->dthetax/ARCSEC,
	  obs->thetay/ARCSEC, obs->dthetay/ARCSEC, obs->obscode);
  return;
}

/* Read an orbit in alpha, beta, gamma format & covar matrix from file.
   Also reads the coordinate frame specifications.
*/
int
read_abg(char *fname, 
	 PBASIS *p, 
	 double **covar)
{
  FILE *fptr;
  int  i;
  char	inbuff[256];

  if (fname==NULL)
    fptr = stdin;
  else if ( (fptr=fopen(fname,"r"))==NULL) {
    fprintf(stderr,"Error opening a/b/g orbit file %s\n",fname);
    return(1);
  }

  /* Skipping comments, read the a/b/g specs*/
  if (fgets_nocomment(inbuff,255,fptr,NULL)==NULL) {
    fprintf(stderr,"Data missing from a/b/g/ data file.\n");
    return(1);
  }

  if (sscanf(inbuff, "%lf %lf %lf %lf %lf %lf",
	     &(p->a),&(p->adot),&(p->b),&(p->bdot),
	     &(p->g),&(p->gdot)) != 6) {
    fprintf(stderr,"Error reading a/b/g data\n");
    return(1);
  }

  for (i=1; i<=6; i++) {
    if (fgets_nocomment(inbuff,255,fptr,NULL)==NULL) {
      fprintf(stderr,"Data missing from a/b/g/ covariance.\n");
      return(1);
    }

    if (sscanf(inbuff, "%lf %lf %lf %lf %lf %lf",
	       &covar[i][1],&covar[i][2],&covar[i][3],
	       &covar[i][4],&covar[i][5],&covar[i][6]) != 6) {
      fprintf(stderr,"Error reading a/b/g covariance\n");
      return(1);
    }
  }

  /* Now read the coordinate system info */
  if (fgets_nocomment(inbuff,255,fptr,NULL)==NULL) {
    fprintf(stderr,"Data missing from a/b/g/ data file.\n");
    return(1);
  }

  if (sscanf(inbuff, "%lf %lf %lf %lf %lf %lf",
	     &lat0, &lon0, &xBary, &yBary, &zBary, &jd0) != 6) {
    fprintf(stderr,"Error reading coord system info\n");
    return(1);
  }
  lat0 *= DTOR;
  lon0 *= DTOR;

  if (fname!=NULL) fclose(fptr);
  return(0);
}

/* Take a set of observations and make a preliminary fit using the
 * linear model of the orbit.  Then fill in zero for gdot, and return
 * the fit and an uncertainty matrix.  The gdot term of uncertainty
 * matrix is set to a nominal value.
 * Note covar is assumed to be 6x6 1-indexed matrix a la Numerical Recipes.
 */
void
prelim_fit(OBSERVATION obsarray[],
	   int nobs,
	   PBASIS *pout,
	   double **covar)
{
  double *beta, *soln, **alpha;
  int	i,j,k;
  double x,y,wtx,wty,*dx,*dy;
  

  beta=dvector(1,6);
  alpha=dmatrix(1,6,1,6);
  soln=dvector(1,6);
  dx=dvector(1,6);
  dy=dvector(1,6);

  /* clear all vectors/matrices*/
  for (i=1; i<=6; i++) {
    beta[i]=soln[i]=0.;
    for (j=1; j<=6; j++) alpha[i][j]=0.;
  }

  /*Collect the requisite sums*/
  for (i=0; i<nobs; i++) {
    wtx = 1./obsarray[i].dthetax;
    wtx *= wtx;
    wty = 1./obsarray[i].dthetay;
    wty *= wty;


    kbo2d_linear(pout,&(obsarray[i]),&x,dx,&y,dy);

    /* Note that the dx[6] and dy[6] terms will only make
     * even the least sense if the pout->g and adot were set to
     * some sensible value beforehand.
     */

    for (j=1; j<=6; j++) {
      beta[j] += obsarray[i].thetax * dx[j] * wtx;
      beta[j] += obsarray[i].thetay * dy[j] * wty;
      for (k=1; k<=j; k++) {
        alpha[j][k] += dx[j]*dx[k]*wtx;
        alpha[j][k] += dy[j]*dy[k]*wty; 
      }
    }
  }

  /* Symmetrize and invert the alpha matrix to give covar.  Note
   * that I'm only going to bother with the first 5 params.
   */
  for (i=1; i<=5; i++)
    for (j=1; j<i; j++)
      alpha[j][i]=alpha[i][j];


  if (invert_matrix(alpha,covar,5)) {
    /* some failure in the inversion...*/
    fprintf(stderr,"Error inverting the alpha matrix\n");
    exit(1);
  }

  /* Now multiply matrices to get the solution vector */
  for (i=1; i<=5; i++) {
    soln[i]=0.;
    for (j=1; j<=5; j++)
      soln[i] += covar[i][j]*beta[j];
  }
   
 /* fill in the PBASIS structure */
  pout->a =    soln[1];
  pout->adot = soln[2];
  pout->b =    soln[3];
  pout->bdot = soln[4];
  pout->g    = soln[5];
  pout->gdot = 0.;              /*assuming that prelim fit has no info here*/

  /*Set the gdot parts of the covariance matrix to nominal values */
  for (i=1; i<6; i++)
    covar[i][6]=covar[6][i]=0.;
  covar[6][6]=0.1*TPI*TPI*pow(pout->g,3.);

  free_dvector(dx,1,6);
  free_dvector(dy,1,6);
  free_dvector(soln,1,6);
  free_dvector(beta,1,6);
  free_dmatrix(alpha,1,6,1,6);
  return;
}

/* Routine to predict position and uncertainty at any time, given
 * a PBASIS fit and a sigma matrix.
 */
void
predict_posn(PBASIS *pin,
             double **covar,
             OBSERVATION *obs,
             double **sigxy)    /*this holds xy error matrix*/
{
  int   i,j,t;
  double *dx,*dy;

  dx=dvector(1,6);
  dy=dvector(1,6);

  /*using input time & obscode, put xy position into OBSERVATION struct*/
  kbo2d(pin,obs,
	&(obs->thetax),dx,&(obs->thetay),dy);

  /* project the covariance matrix */
  /* skip if not desired */
  if (sigxy!=NULL && covar!=NULL) {
    sigxy[1][1]=sigxy[1][2]=sigxy[2][2]=0;
    for (i=1; i<=6; i++) { 
      for (j=1; j<=6; j++) {
	sigxy[1][1] += dx[i]*covar[i][j]*dx[j];
	sigxy[1][2] += dx[i]*covar[i][j]*dy[j];
	sigxy[2][2] += dy[i]*covar[i][j]*dy[j];
      }
    }
    sigxy[2][1]=sigxy[1][2];
  }

  free_dvector(dx,1,6);
  free_dvector(dy,1,6);
  return;
}

/* Invert a double matrix, 1-indexed of size dim
 * from Numerical Recipes.  Input matrix is destroyed.
 */
int
invert_matrix(double **in, double **out,
              int dim)
{
  extern void ludcmp(double **a, int n, int *indx, double *d);
  extern void ludcmp(double **a, int n, int *indx, double *d);

  int   *indx,i,j;
  double *tvec,det;

  tvec = dvector(1,dim);
  indx = ivector(1,dim);
  ludcmp(in,dim,indx,&det);

  for (j=1; j<=dim; j++) {
    for (i=1; i<=dim; i++) tvec[i]=0.;
    tvec[j] = 1.0;
    lubksb(in,dim,indx,tvec);
    for (i=1; i<=dim; i++) out[i][j]=tvec[i];
  }

  free_ivector(indx,1,6);
  free_dvector(tvec,1,6);
  return(0);
}

void
print_matrix(FILE *fptr, double **matrix, int xdim, int ydim)
{
  int   i,j;
  for (i=1; i<=ydim; i++) {
    for (j=1; j<=xdim; j++) {
      fprintf(fptr,"%11.4e ",matrix[i][j]);
    }
    fprintf(fptr,"\n");
  }
  return;
}

/* Function to return xyz coords of observatory in the current
 * standard coordinate system.
 */
void
earth3d(double t,	/* time is in years here */
	int obscode,
	double *x, double *y, double *z)
{
  double x1,y1,z1,xTelEq,yTelEq,zTelEq;
  double xec, yec, zec; 

  /* get observatory posn wrt barycenter */
  earth_ssbary(t/DAY+jd0, obscode, &x1, &y1, &z1);
  /* convert to tangent-point coord system */
  /* via ecliptic */
  xyz_eq_to_ec(x1, y1, z1, &xec, &yec, &zec,NULL);
  xyz_ec_to_proj(xec, yec, zec, x, y, z, lat0, lon0, NULL);

  /* Translate to our origin */
  *x += xBary;
  *y += yBary;
  *z += zBary;

  return;
}


/* Routine to give the t=0 position & velocity of KBO in
 * ecliptic system centered at SSBARY */
void
pbasis_to_bary(PBASIS *p,
	       XVBASIS *xv,
	       double **partials)
{
  double xProj,yProj,zProj,xdotProj,ydotProj,zdotProj; 
  double z0; 
  double **dProj_dp, **dEc_dProj;	/*intermediate partial matrices */
  int  i,j;

  if (partials!=NULL) {
    dProj_dp = dmatrix(1,6,1,6);
    dEc_dProj = dmatrix(1,6,1,6);
    for (i=1; i<=6; i++)
      for (j=1; j<=6; j++)
	dProj_dp[i][j]=dEc_dProj[i][j]=0.;
  }

  z0 = 1./p->g;
  xProj=p->a*z0;
  yProj=p->b*z0;
  zProj=z0;
  xdotProj=p->adot *z0;
  ydotProj=p->bdot *z0;
  zdotProj=p->gdot *z0;

  if (partials!=NULL) {
    dProj_dp[1][1] = z0;
    dProj_dp[1][5] = -p->a*z0*z0;
    dProj_dp[2][3] = z0;
    dProj_dp[2][5] = -p->b*z0*z0;
    dProj_dp[3][5] = -z0*z0;
    dProj_dp[4][2] = z0;
    dProj_dp[4][5] = -p->adot*z0*z0;
    dProj_dp[5][4] = z0;
    dProj_dp[5][5] = -p->bdot*z0*z0;
    dProj_dp[6][6] = z0;
    dProj_dp[6][5] = -p->gdot*z0*z0;
  }

  /* Bring position into ecliptic system relative to SSBARY */
  if (partials!=NULL)
    xyz_proj_to_ec(xProj-xBary, yProj-yBary, zProj-zBary,
		   &(xv->x), &(xv->y), &(xv->z),
		   lat0, lon0, dEc_dProj);
  else 
    xyz_proj_to_ec(xProj-xBary, yProj-yBary, zProj-zBary,
		   &(xv->x), &(xv->y), &(xv->z),
		   lat0, lon0, NULL);

  /* Repeat for velocity */
  xyz_proj_to_ec(xdotProj, ydotProj, zdotProj,
		 &(xv->xdot), &(xv->ydot), &(xv->zdot),
		 lat0, lon0, NULL);

  /* Calculate the overall Jacobian if desired*/
  if (partials!=NULL) {
    /* Replicate the spatial part of dEc_dProj into velocities */
    for (i=4; i<=6; i++)
      for (j=4; j<=6; j++)
	dEc_dProj[i][j] = dEc_dProj[i-3][j-3];

    matrix_multiply(dEc_dProj,dProj_dp,partials,6,6,6,6);

    free_dmatrix(dProj_dp,1,6,1,6);
    free_dmatrix(dEc_dProj,1,6,1,6);
  }

  return;   
}

/* Multiply matrix m1 x m2.  m1 is x1 cols by y1 rows, etc.
   A dumb implementation.  Note that SECOND index of
   matrices are varying most rapidly, along rows.
*/
void
matrix_multiply(double **m1, double **m2,
		double **mout,
		int x1, int y1, int x2, int y2)
{
  int	i,j,k;
  double	sum;

  if (x1!=y2) {
    fprintf(stderr,"Trying to multiply mismatched matrices, "
	    " %d x %d  times %d x %d\n", y1,x1,y2,x2);
    exit(1);
  }
  for (i=1; i<=y1; i++)
    for (k=1; k<=x2; k++) {
      sum = 0.;
      for (j=1; j<=y2; j++)
	sum += m1[i][j]*m2[j][k];
      mout[i][k] = sum;
    }

  return;
}

/* Remap covariance matrix from one basis to another, using
   the partial deriv matrix **derivs.
   kin = dimension on input, kout = dimension on output.
   Just calculates deriv_T x covar_in x deriv.
*/
void
covar_map(double **covar_in, 
	  double **derivs, 
	  double **covar_out,
	  int kin, int kout)
{
  int	i,j,m,n;
  double	sum;

  for (i=1; i<=kout; i++)
    for (j=1; j<=kout; j++) {
      sum = 0.;
      for (m=1; m<=kin; m++)
	for (n=1; n<=kin; n++)
	sum += derivs[i][m]*derivs[j][n]*covar_in[m][n];
      covar_out[i][j] = sum;
    }

  return;
}
 
/* calculate acceleration given barycentric coords x */
/* This version just does the barycentric approximation. */
void
accel(double *xx,
      double *a,
      double t)
{
  double x,y,z,acc,r2;
  x = xx[0] - xBary;
  y = xx[1] - yBary;
  z = xx[2] - zBary;
  r2 = x*x+y*y+z*z;
  acc = -GM*SSMASS*pow( r2 , -1.5);
  a[0] = x*acc;
  a[1] = y*acc;
  a[2] = z*acc;
  return;
}

/* calculate acceleration given barycentric coords x */
/* Here is a more sophisticated version that includes perturbations*/
/**** ??? Note that we could cache the planet positions if this turns
*** out to be slower than desired.
***/
void
accel_nbody(double *xx,
	    double *a,
	    double t)
{
  double x,y,z,acc,r2;
  double xb, yb, zb;
  static int bodies[]={SUN, JUPITER, SATURN, URANUS, NEPTUNE};
  static int nbodies = 5; 
  static double mass[]={1.000006, 9.548e-4, 2.859e-4, 4.37e-5, 5.15e-5};
  /* Note Sun mass includes the inner planets*/
  int i,j;

  for (j=0; j<3; j++) a[j]=0.;

  for (i=0; i<nbodies; i++) {
    body3d(t, bodies[i], &xb, &yb, &zb, NULL);
    x = xx[0] - xb;
    y = xx[1] - yb;
    z = xx[2] - zb;
    r2 = x*x+y*y+z*z;
    acc = mass[i]*pow( r2 , -1.5);
    a[0] += x*acc;
    a[1] += y*acc;
    a[2] += z*acc;
  }

  for (j=0; j<3; j++) a[j]*= -GM;

  return;
}

/* Give the KBO's 3-d position, along with derivatives. */
/* Here is a leapfrog-integrator version.*/
void    
kbo3d(PBASIS *pin,
      double t,
      double *xout,
      double *vout,
      double dx[],
      double dy[],
      double dz[])
{
  int   i, restart=1;
  static double x[3],v[3],a[3], tx, tv, z0;
  static PBASIS psave;
  static int init=0;
  static int tdir;
  double t1,dt,dtv;
  double tstep=20.*DAY; 

  /* decide whether we need to reset integrator to t=0*/
  if (!init) {
    init = 1;
    restart = 1;
  } else if (pin->a==psave.a && pin->adot==psave.adot &&
	     pin->b==psave.b && pin->bdot==psave.bdot &&
	     pin->g==psave.g && pin->gdot==psave.gdot) {
    restart = 0;
  } else {
    /* mismatch in parameters */
    restart=1;
  }

  /* also restart if we'd need to reverse the integration */
  if (restart==0) {
    if (tdir>0 && t<tx-0.6*TSTEP) restart=1;
    if (tdir<0 && t>tx+0.6*TSTEP) restart=1;
  }

  if (restart) {
    /* initial time & velocity */
    z0 = 1./ pin->g;
    x[0] = pin->a *z0;
    x[1] = pin->b *z0;
    x[2] = z0;
    v[0] = pin->adot*z0;
    v[1] = pin->bdot*z0;
    v[2] = pin->gdot*z0;
    tx = tv = 0.;
    psave.a = pin->a;
    psave.adot = pin->adot;
    psave.b = pin->b;
    psave.bdot = pin->bdot;
    psave.g = pin->g;
    psave.gdot = pin->gdot;

    /* choose direction of integrator */
    if (t>=0.) tdir=+1;
    else tdir=-1;
  }

  /* leap-frog until tx is within TSTEP/2 of target time*/
  for ( ; tdir*(t-tx)>TSTEP/2.; tx = t1) {
    t1 = tx + tdir*TSTEP;
    dt = t1 - tx;
    
    /* jump velocity to middle of time step */
    accel_nbody(x,a,tx);
    /*accel(x,a,tx);*/
    dtv = 0.5*dt + (tx-tv);
    for (i=0; i<3; i++) v[i] += dtv*a[i];
    tv += dtv;

    /* Jump position using this velocity */
    for (i=0; i<3; i++) x[i] += dt*v[i];
  }

  /* Now take the last step from tx to t */
  accel_nbody(x,a,tx);
  /*accel(x,a,tx);*/
  for (i=0; i<3; i++) {
    xout[i] = x[i] + v[i]*(t-tx) + a[i]*(t-tx)*(0.5*(t+tx)-tv);
    vout[i] = v[i] + a[i]*(t-tv);
  }
  /* x and y derivatives w.r.t. parameters  - inertial approx ok here*/
  if (dx== NULL || dy==NULL || dx==NULL) return;

  for (i=1; i<=6; i++) dx[i]=dy[i]=dz[i]=0.;
  dx[1] = dy[3] = z0;
  dx[2] = dy[4] = t * z0;
  dx[5] = xout[0] * (-z0);
  dy[5] = xout[1] * (-z0);
  dz[5] = xout[2] * (-z0);
  dz[6] = t*z0;

  return;
}

/* Function to return xyz coords of a JPL ephemeris body in
 * standard coordinate system.
 */
void
body3d(double t,	/* time is in years here */
       int body,
       double *x, double *y, double *z,
       double *vxyz)
{
  double xxx[3];
  double xec, yec, zec; 

  /* get observatory posn wrt barycenter */
  bodycenter_ssbary(t/DAY+jd0, xxx, body, vxyz);
  /*fprintf(stderr,"bodycenter for %d: %f %f %f\n",body,xxx[0],xxx[1],xxx[2]);*/
  /* convert to tangent-point coord system */
  /* via ecliptic */
  xyz_eq_to_ec(xxx[0], xxx[1], xxx[2], &xec, &yec, &zec,NULL);
  xyz_ec_to_proj(xec, yec, zec, x, y, z, lat0, lon0, NULL);

  /* Translate to our origin */
  *x += xBary;
  *y += yBary;
  *z += zBary;

  /* Rotate velocity if it has bee requested: */
  if (vxyz != NULL) {
    xyz_eq_to_ec(vxyz[0], vxyz[1], vxyz[2], &xec, &yec, &zec,NULL);
    xyz_ec_to_proj(xec, yec, zec, vxyz, vxyz+1, vxyz+2, lat0, lon0, NULL);
  }
    
  return;
}

/* Return the elongation of an observation re barycenter */
double
elongation(OBSERVATION *obs)
{
  double cosb,ee;

  if (obs->xe < -9.) {
    earth3d(obs->obstime, obs->obscode, &(obs->xe),&(obs->ye),&(obs->ze));
  }
  /* dot the observation direction into the (Earth->bary) vector*/
  cosb = obs->thetax*(xBary-obs->xe) + obs->thetay*(yBary-obs->ye) +
    (zBary - obs->ze);
  cosb /= sqrt(1+obs->thetax*obs->thetax+obs->thetay*obs->thetay);
  ee = (xBary-obs->xe)*(xBary-obs->xe);
  ee += (yBary-obs->ye)*(yBary-obs->ye);
  ee += (zBary-obs->ze)*(zBary-obs->ze);
  cosb /= sqrt(ee);
  return acos(cosb);
}

/* Return the difference in ecliptic longitude between object & bary
 * (e.g., 180 degrees is opposition).
 * gives Sun - Target so it increases about 1 deg/day.
 */
double
opposition_angle(OBSERVATION *obs)
{
  double lat, lonkbo, lonbary;
  double xx, yy, zz, xec, yec, zec;

  if (obs->xe < -9.) {
    earth3d(obs->obstime, obs->obscode, &(obs->xe),&(obs->ye),&(obs->ze));
  }

  /* Get barycenter ecliptic longitude */
  xx = xBary - obs->xe;
  yy = yBary - obs->ye;
  zz = zBary - obs->ze;
  xyz_proj_to_ec(xx, yy, zz, &xec, &yec, &zec, lat0, lon0, NULL);
  lonbary = atan2(yec,xec);

  /* And the observation ecliptic longitude*/
  proj_to_ec(obs->thetax, obs->thetay, &lat, &lonkbo, lat0, lon0, NULL);

  lonkbo = lonbary - lonkbo;
  while (lonkbo<0) lonkbo += 2*PI;

  return lonkbo;
}

/* Return the angle btwn zenith (anti-earth) direction and target*/
double
zenith_angle(OBSERVATION *obs)
{
  double xobs, yobs, zobs, r;
  double xec, yec, zec, cosb;

  /* Get the anti-Earth vector (in ICRS) for this observation */
  observatory_geocenter(obs->obstime/DAY+jd0, obs->obscode, 
			&xobs, &yobs, &zobs);
  r = sqrt(xobs*xobs+yobs*yobs+zobs*zobs);
  if (r<=0.) {
    fprintf(stderr,"Non-positive geocentric radius in zenith_angle()\n");
    exit(1);
  }
  xobs /=r; yobs/=r; zobs/=r;
  /* Rotate this ICRS vector into ecliptic, then projected coords */
  xyz_eq_to_ec(xobs, yobs, zobs, &xec, &yec, &zec,NULL);
  xyz_ec_to_proj(xec, yec, zec, &xobs, &yobs, &zobs, lat0, lon0, NULL);

  /* dot the observation direction into the (Earth->bary) vector*/
  cosb = obs->thetax*xobs + obs->thetay*yobs + zobs;
  cosb /= sqrt(1+obs->thetax*obs->thetax+obs->thetay*obs->thetay);
  return acos(cosb);
}

/* return true (1) if target is visible (above horizon) */
int
is_visible(OBSERVATION *obs)
{
  return zenith_angle(obs) < zenith_horizon(obs->obscode);
}

#include <sys/time.h>

void
fake_observation(PBASIS *p, 
		 OBSERVATION *obs)
{
  static long idum=-1;
  float gasdev(long *idum);

  /* seed the random number generator*/
  if (idum<0) {
    struct timeval tp;
    gettimeofday(&tp,NULL);
    idum = -tp.tv_usec;
  }

  kbo2d(p,obs,&(obs->thetax),NULL,&(obs->thetay),NULL);

  /* Add measurement noise to the positions */
  obs->thetax += obs->dthetax*gasdev(&idum);
  obs->thetay += obs->dthetay*gasdev(&idum);

  return;
}

/* Parse command-line looking for the standard options.
 * iarg is first arg to examine, and is returned as the
 * first non-option argument.
 * Returns 1 on parse error.
 * currently understands m, j, o, and v options.
 */
#define BUFFSIZE 512
int
read_options(int* iarg, 
	     int argc,
	     char *argv[])
{
  double d;
  for ( ; *iarg<argc; (*iarg)++) {
    if (argv[*iarg][0]!='-')
      return(0);	/*quit, now at non-option argument.*/

    if (strcasecmp(argv[*iarg]+1,"m")==0) {
      /* MPC error size */
      d = atof(argv[++(*iarg)]);
      if (d<=0.) {
	fprintf(stderr,"Bad MPC error spec %s\n",argv[*iarg]);
	return(1);
      }
      set_mpc_dtheta(d);

    } else if (strcasecmp(argv[*iarg]+1,"o")==0) {
      /* observatory file name*/
      set_observatory_file(argv[++(*iarg)]);

    } else if (strcasecmp(argv[*iarg]+1,"j")==0) {
      /* ephemeris file name*/
      set_ephem_file(argv[++(*iarg)]);

    } else if (strcasecmp(argv[*iarg]+1,"v")==0) {
      /* version number request - quit after */
      printf("KBO orbit-fitting software Release is $Name: Release2_0 $\n");
      printf("Contact Gary Bernstein for information on use,\n");
      printf("  problems, and credit for this software.\n");
      exit(1);

    } else {
      fprintf(stderr,"Unknown option %s\n",argv[*iarg]);
      return(1);
    }
  }
  return(0);
}
