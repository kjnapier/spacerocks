/* 	$Id: orbfit.h,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */
/* Definitions for the orbit-fitting software */
#include "nrutil.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
 
/*maximum number of observations*/
#define MAXOBS  4000

/*Default name of binary DE423 ephemeris data:*/
#define DEFAULT_EPHEM_FILE "binEphem.423" 
/*Or override the above by looking for filename under
 * this environment variable:
 */
#define EPHEM_ENVIRON "ORBIT_EPHEMERIS"

/*Name of observatory location file:*/
#define DEFAULT_OBSERVATORY_FILE "observatories.dat"
/*Or override the above by looking for filename under
 * this environment variable:
 */
#define OBS_ENVIRON "ORBIT_OBSERVATORIES"

/*maximum number of observatories in data file*/
#define MAX_SITES  400    

/* Some magic numbers for observatories */
#define OBSCODE_SSBARY 000
#define OBSCODE_GEOCENTER 500
/* obscodes above following are taken to be orbiting earth,
 * below this value they are affixed to earth */
#define OBSCODE_ORBITAL 2000

/*uncertainty assumed for MPC-format observations (arcsec)*/
#define DEFAULT_DTHETA 0.3

#ifndef PI
#define PI      3.14159265359
#endif
#define TPI     (2.*PI)
#define DTOR	(PI/180.)
#define GM      (4.*PI*PI/1.0000378)    /*solar gravitation*/
#define SSMASS  1.00134       /*Total SS mass*/
#define ARCSEC	(PI/180./3600.)
#define DAY	(1./365.25)	/*Julian day, 86400 s*/
#define SPEED_OF_LIGHT  63241.06515	/*in AU/YR*/
#define ECL	(23.43928*PI/180.)	/*Obliquity of ecliptic at J2000*/

#define TSTEP	(20.*DAY)	/*Time step for orbit integrator*/
#define EARTHMASS 3.00349e-6	/*Earth mass relative to Sun (for 
				 *orbiting telescopes) */

/* an orbit as specified by its posn & velocity in fitting params */
typedef struct {
  double a;             /*alpha, angular position at t=0*/
  double adot;          /* (dx/dt) / z at t=0 */
  double b,bdot;        /* same 2 for y direction */
  double g,gdot;        /* 1/z and (dz/dt)/z of KBO at t=0 */
} PBASIS;

/* an orbit as specified by actual spatial position and velocity: */
typedef struct {
  double x, y, z;       /*3-space position */
  double xdot, ydot, zdot;       /*3-space velocity */
  double jd0;		/* Time at this phase space location*/
} XVBASIS;
 
/* an orbit as specified by the usual orbital parameters */
typedef struct {
  double a,e,i;                 /*semi-major axis, ellipticity, inclination*/
  double lan, aop, T;       /*long. of ascending node, arg of perihelion, and
                                  time from periapse of object */
  double ma;                 /* mean anomaly */
} ORBIT;
 
/* From skycalc, altered a bit: */
struct date_time
   {
	int y;
	int mo;
	float d;
	float h;
	float mn;
	float s;
   };

/* Data of an individual observation */
/* Notice that in practice we will have to define an axis for the
 * coordinate system, origin of coordinate system, and time zeropoint.
 * These will for convenience be taken to correspond to first observation
 * in most cases.
 */
typedef struct {
  double        thetax,dthetax; /*x-direction (ecliptic) position
                                  and uncertainty*/
  double        thetay,dthetay; /*y-direction*/
  double        obstime;        /*time of observation (assumed years)*/
  int		obscode;        /*Observatory site ID (from MPC list)*/
  double	xe, ye, ze;	/*Position of observatory at obstime */
} OBSERVATION;
 
/* Declare the variables which define the coordinate system */
extern double	lat0, lon0;	/* ecliptic lat & lon of tangent point */
extern double	xBary, yBary, zBary;	/*Posn of barycenter in our system*/
extern double	jd0;		/* Zeropoint of time scale */

/* Some functions that we'll need */
 
/* give 3-space coords of KBO and derivs w.r.t. parameters */
extern  void    
kbo3d(PBASIS *pin,
      double t,
      double *xout,
      double *yout,
      double dx[],
      double dy[],
      double dz[]);

/* 3-space position of Earth */
void
earth3d(double t,	/* time is in years here */
	int obscode,
	double *x, double *y, double *z);
/*ICRS vector from SSBary to observatory location */
void
earth_ssbary(double	jd,
	     int	obscode,
	     double *x, double *y, double *z);

void
body3d(double t,	/* time is in years here */
       int body,
       double *x, double *y, double *z,
       double *vxyz);

/*ICRS vector from SSBary to center of some body*/
void
bodycenter_ssbary(double jd,
		  double *xyz,
		  int body,
		  double *vxyz);

/*ICRS vector from SSBary to geocenter */
void
geocenter_ssbary(double jd,
		 double *xyz);

/* ICRS vector from geocenter to observatory location*/
void
observatory_geocenter(double jd,
		      int obscode,
		      double *xobs,
		      double *yobs,
		      double *zobs);

/* Read the observatory location from fname.  fname=NULL looks
 * in default/environment places.*/
void
read_observatories(char *fname);

/* Generic help dumper */
void
print_help(char *h[]);


/* Angle (radians) from zenith to the limb/horizon */
double
zenith_horizon(int obscode);

/* Angle (radians) from zenith to observation direction */
double
zenith_angle(OBSERVATION *o);

/* return true if target is observable (above horizon) */
int
is_visible(OBSERVATION *o);

/* tangent-plane angular position of KBO and derivs */
extern  void
kbo2d(PBASIS *pin, 
      OBSERVATION *obs,
      double *x, double dx[],
      double *y, double dy[]);
 
/* linearized version of the 2d position, ignores gdot term*/
extern  void
kbo2d_linear(PBASIS *pin,
	     OBSERVATION *obs,
	     double *x, double dx[],
             double *y, double dy[]);
 
/* map from PBASIS to ORBIT format, and derivative matrix*/
/*
 *extern int    abg_to_aei(PBASIS *pin, ORBIT *orbout,
 *                          double **pderivs);
 */
/* extract observations from a string */
int
scan_observation(char *inbuff,
		 OBSERVATION *obs);
/* read RA/DEC from file & set up coord systems */
int
read_radec(OBSERVATION obsarray[], 
	   char *fname, 
	   int *nobs);
/* Rolls all the i/o and fitting together */
void
fit_radec(char *fname, int *nobservations, double *chisqfit, int *ndof, PBASIS *pbasis, ORBIT *orb);
/* reset astrometric error assigned to MPC obs.*/
void
set_mpc_dtheta(double d);

/* set filename for ephemeris or observatory data */
void
set_ephem_file(char *fname);
void
set_observatory_file(char *fname);

/* Read these options */
int
read_options(int* iarg, 
	     int argc,
	     char *argv[]);

/* preliminary fit to observations */
void
prelim_fit(OBSERVATION obsarray[],
	   int nobs,
	   PBASIS *pout,
	   double **covar);

/* Routine to predict position and uncertainty at any time, given
 * a PBASIS fit and a sigma matrix.
 */
void
predict_posn(PBASIS *pin,
             double **covar,
             OBSERVATION *obs,
             double **sigxy);

/* Take a PBASIS orbit and uncertainty and convert it to a traditional
 * ORBIT basis and uncertainty.
 */

/* void
 * predict_aei(PBASIS *pin, double **covarp,
 *	    ORBIT  *orbout, double **covarq);
 */

void
print_matrix(FILE *fptr, double **matrix, int xdim, int ydim);

char *
fgets_nocomment(char *buffer, int length, FILE *fptr, FILE *fpout);
/* Coordinate transformation routines: */
void
eq_to_ec( double ra_eq,   double dec_eq,
	  double *lat_ec, double *lon_ec,
	  double **partials);
void
xyz_eq_to_ec(double x_eq, double y_eq, double z_eq,
	     double *x_ec, double *y_ec, double *z_ec,
	     double **partials);
void
ec_to_eq( double lat_ec,  double lon_ec,
	  double *ra_eq,  double *dec_eq,
	  double **partials);
void
xyz_ec_to_eq(double x_ec, double y_ec, double z_ec,
	     double *x_eq, double *y_eq, double *z_eq,
	     double **partials);
void
ec_to_proj(double lat_ec,   double lon_ec,
	   double *x_proj,  double *y_proj,
	   double lat0,     double lon0,
	   double **partials);
void
proj_to_ec(double x_proj,   double y_proj,
	   double *lat_ec,  double *lon_ec,
	   double lat0,	    double lon0,
	   double **partials);
void
xyz_ec_to_proj( double x_ec, double y_ec, double z_ec,
		double *x_p, double *y_p, double *z_p,
		double lat0, double lon0,
		double **partials);
void
xyz_proj_to_ec( double x_p, double y_p, double z_p,
		double *x_ec, double *y_ec, double *z_ec,
		double lat0, double lon0,
		double **partials);
void
pbasis_to_bary(PBASIS *p,
	       XVBASIS *xv,
	       double **partials);
void
matrix_multiply(double **m1, double **m2,
		double **mout,
		int x1, int y1, int x2, int y2);
void
orbitElements(XVBASIS *xv,
	      ORBIT  *orb);
void
elements_to_xv(ORBIT *o,
	       double jd,
	       XVBASIS *xv);
void
elements_to_pbasis(ORBIT *o,
		   double jd,
		   int obscode,
		   PBASIS *p);
void
covar_map(double **covar_in, 
	  double **derivs, 
	  double **covar_out,
	  int kin, int kout);
double 
date_to_jd(struct date_time date);
void
aei_derivs( XVBASIS *xv,
	    double **daei_dxv);

double
elongation(OBSERVATION *obs);
double
opposition_angle(OBSERVATION *obs);

void
fake_observation(PBASIS *p, 
		 OBSERVATION *obs);

void
flatten_cov(double **cov, int ndim, double *cov1d);

void
unflatten_cov(double *cov1d, int ndim, double **cov);

void
add_to_obsarray(OBSERVATION obsarray[], int iobs, OBSERVATION obs);

