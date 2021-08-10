/* 	$Id: ephem_earth.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */
#ifndef lint
static char vcid[] = "$Id: ephem_earth.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/*** ephem_earth.c  - I've changed the below to include a routine explicitly
*** returning the location of earth geocenter relative to SSBARY.  Also
*** have eliminated the nutation & libration routines.
*** 8/9/99 gmb
*** add 64s difference between UT (usually used in position tables) and TDB
*** (used in this ephemeris).
***/
#include <string.h>
#include "orbfit.h"
#define TDBOFFSET (64./86400.)
#define FNAMESIZE 256
#define BUFFSIZE  512
/* structures used to keep observatory information: */
static int nsites=0;
typedef struct {
  int code;
  double	lon;		/*observatory longitude in hours*/
  double	lat;		/*observatory latitude in degrees*/
  double	altitude;	/*observatory altitude in meters */
  char	name[80];
} SITE;
static SITE sitelist[MAX_SITES];

/* And likewise some information for orbiting observatories: */
static int nspacecraft=0;
typedef struct {
  int code;
  double	i;		/*inclination of orbit (degrees)*/
  double	P;		/*sidereal orbital period, (days)*/
  double	precess;	/*period for orbital precession (days)*/
  double	jd0;		/*time of zero orbit phase*/
  double	ra0;		/*RA of orbit pole at jd0 (degrees)*/
  double	a;		/*semi-major axis of orbit (m)*/
  char	name[80];
} SPACECRAFT;
static SPACECRAFT spacecraftlist[MAX_SITES];

char  observatory_file[FNAMESIZE]="";
void
set_observatory_file(char *fname) {
  strncpy(observatory_file, fname, FNAMESIZE-1);
  observatory_file[FNAMESIZE-1]=0;
}

char  ephem_file[FNAMESIZE]="";
void
set_ephem_file(char *fname) {
  strncpy(ephem_file, fname, FNAMESIZE-1);
  ephem_file[FNAMESIZE-1]=0;
}



/******************************************************************************/
/**                                                                          **/
/**  SOURCE FILE: ephem_read.c                                               **/
/**                                                                          **/
/**     This file contains a set of functions and global variables that      **/
/**     implements an ephemeris server program. Client programs can use      **/
/**     use this server to access ephemeris data by calling one of the       **/
/**     following functions:                                                 **/
/**                                                                          **/
/**        Interpolate_Libration -- returns lunar libration angles           **/
/**        Interpolate_Nutation  -- returns (terrestrial) nutation angles    **/
/**        Interpolate_Position  -- returns position of planet               **/
/**        Interpolate_State     -- returns position and velocity of planet  **/
/**                                                                          **/
/**     Note that client programs must make one, and only one, call to the   **/
/**     function Initialize_Ephemeris before any of these other functions    **/ 
/**     be used. After this call any of the above functions can be called    **/
/**     in any sequence.                                                     **/
/**                                                                          **/
/**  Programmer: David Hoffman/EG5                                           **/
/**              NASA, Johnson Space Center                                  **/
/**              Houston, TX 77058                                           **/
/**              e-mail: david.a.hoffman1@jsc.nasa.gov                       **/
/**                                                                          **/
/******************************************************************************/

#include <stdio.h>
#include <math.h>
#ifndef TYPES_DEFINED
#include "ephem_types.h"
#endif

#ifndef PI
#define PI 3.1415926535
#endif

/**==========================================================================**/
/**  Global Variables                                                        **/
/**==========================================================================**/

   static headOneType  H1;
   static headTwoType  H2;
   static recOneType   R1;
   static FILE        *Ephemeris_File;
   static double       Coeff_Array[ARRAY_SIZE] , T_beg , T_end , T_span;

   static int Debug = FALSE;             /* Generates detailed output if true */

/**==========================================================================**/
/**  Read_Coefficients                                                       **/
/**                                                                          **/
/**     This function is used by the functions below to read an array of     **/
/**     Tchebeychev coefficients from a binary ephemeris data file.          **/
/**                                                                          **/
/**  Input: Desired record time.                                             **/
/**                                                                          **/
/**  Output: None.                                                           **/
/**                                                                          **/
/**==========================================================================**/

void Read_Coefficients( double Time )
{
  double  T_delta = 0.0;
  long     Offset  =  0 ;		/*** ??? change to long 8/9/99 ***/

  /*--------------------------------------------------------------------------*/
  /*  Find ephemeris data that record contains input time. Note that one, and */
  /*  only one, of the following conditional statements will be true (if both */
  /*  were false, this function would not have been called).                  */
  /*--------------------------------------------------------------------------*/

  if ( Time < T_beg )                    /* Compute backwards location offset */
     {
       T_delta = T_beg - Time;
       Offset  = (int) -ceil(T_delta/T_span);	/***Needed negative sign here???*/
     }

  if ( Time > T_end )                    /* Compute forewards location offset */
     {
       T_delta = Time - T_end;
       Offset  = (int) ceil(T_delta/T_span);
     }

  /*--------------------------------------------------------------------------*/
  /*  Retrieve ephemeris data from new record.                                */
  /*--------------------------------------------------------------------------*/

  fseek(Ephemeris_File,(Offset-1)*ARRAY_SIZE*sizeof(double),SEEK_CUR);
  fread(&Coeff_Array,sizeof(double),ARRAY_SIZE,Ephemeris_File);
  
  T_beg  = Coeff_Array[0];
  T_end  = Coeff_Array[1];
  T_span = T_end - T_beg;

  if (Time < T_beg || Time > T_end) {
   /* fprintf(stderr,"JD %.2f is out of range of ephemeris file\n",Time); */
 /*   exit(1); */
    Time = T_beg;  /* DG hack to avoid exiting in this situation */
  }
  /*--------------------------------------------------------------------------*/
  /*  Debug print (optional)                                                  */
  /*--------------------------------------------------------------------------*/

  if ( Debug && T_beg < 2400000.0) 
     {
 /*      printf("\n  In: Read_Coefficients \n");
       printf("\n      ARRAY_SIZE = %4d",ARRAY_SIZE);
       printf("\n      Offset  = %3d",Offset);
       printf("\n      T_delta = %7.3f",T_delta);
*/
       printf("\n      T_Beg   = %7.3f",T_beg);
  /*
       printf("\n      T_End   = %7.3f",T_end);
       printf("\n      T_Span  = %7.3f\n\n",T_span);
*/
     }

}

/**==========================================================================**/
/**  Initialize_Ephemeris                                                    **/
/**                                                                          **/
/**     This function must be called once by any program that accesses the   **/
/**     ephemeris data. It opens the ephemeris data file, reads the header   **/
/**     data, loads the first coefficient record into a global array, then   **/
/**     returns a status code that indicates whether or not all of this was  **/
/**     done successfully.                                                   **/
/**                                                                          **/
/**  Input: A character string giving the name of an ephemeris data file.    **/
/**                                                                          **/
/**  Returns: An integer status code.                                        **/
/**                                                                          **/
/**==========================================================================**/
/** This routine no longer takes filename as argument, it goes looking
 ** for the file itself, as environment-specified file or as default in
 ** this directory
 **/
int Initialize_Ephemeris()
{
  int headerID;
  char fileName[FNAMESIZE];
  /*** gmb: don't duplicate this call */
  static init=0;
  if (init) return SUCCESS;
  init = 1;

  /*--------------------------------------------------------------------------*/
  /*  Open ephemeris file.                                                    */
  /*--------------------------------------------------------------------------*/

  /** use previously specified filename,
   ** or environment-specified file, or the default filename, in that order
   */
  if (strlen(ephem_file)>0)
    strncpy(fileName, ephem_file, FNAMESIZE-1);
  else if (getenv(EPHEM_ENVIRON)!=NULL)
    strncpy(fileName, getenv(EPHEM_ENVIRON), FNAMESIZE-1);
  else
    strncpy(fileName, DEFAULT_EPHEM_FILE, FNAMESIZE-1);


  fileName[FNAMESIZE-1]=0;

  
  Ephemeris_File = fopen(fileName,"r");

  /*--------------------------------------------------------------------------*/
  /*  Read header & first coefficient array, then return status code.         */
  /*--------------------------------------------------------------------------*/

  if ( Ephemeris_File == NULL ) /*........................No need to continue */
     {
       printf("\n Unable to open ephemeris file: %s.\n",fileName);
       return FAILURE;
     }
  else 
     { /*.................Read first three header records from ephemeris file */
         
       fread(&H1,sizeof(double),ARRAY_SIZE,Ephemeris_File);
       fread(&H2,sizeof(double),ARRAY_SIZE,Ephemeris_File);
       fread(&Coeff_Array,sizeof(double),ARRAY_SIZE,Ephemeris_File);
       
       /*...............................Store header data in global variables */
       
       R1 = H1.data;
              
       /*..........................................Set current time variables */

       T_beg  = Coeff_Array[0];
       T_end  = Coeff_Array[1];
       T_span = T_end - T_beg;

       /*..............................Convert header ephemeris ID to integer */

       headerID = (int) R1.DENUM;
       
       /*..............................................Debug Print (optional) */

       if ( Debug ) 
          {
            printf("\n  In: Initialize_Ephemeris \n");
            printf("\n      ARRAY_SIZE = %4d",ARRAY_SIZE);
            printf("\n      headerID   = %3d",headerID);
            printf("\n      T_Beg      = %7.3f",T_beg);
            printf("\n      T_End      = %7.3f",T_end);
            printf("\n      T_Span     = %7.3f\n\n",T_span);
          }

       /*..................................................Return status code */
       
       if ( headerID == EPHEMERIS ) 
          {
            return SUCCESS;
          }
       else 
          {
            printf("\n Opened wrong file: %s",fileName);
            printf(" for ephemeris: %d.\n",EPHEMERIS);
            return FAILURE;
          }
     }
}

/**==========================================================================**/
/**  Interpolate_Position                                                    **/
/**                                                                          **/
/**     This function computes a position vector for a selected planetary    **/
/**     body from Chebyshev coefficients read in from an ephemeris data      **/
/**     file. These coefficients are read from the data file by calling      **/
/**     the function Read_Coefficients (when necessary).                     **/
/**                                                                          **/
/**  Inputs:                                                                 **/
/**     Time     -- Time for which position is desired (Julian Date).        **/
/**     Target   -- Solar system body for which position is desired.         **/
/**     Position -- Pointer to external array to receive the position.       **/
/**                                                                          **/
/**  Returns: Nothing explicitly.                                            **/
/**                                                                          **/
/**==========================================================================**/

void Interpolate_Position( double Time , int Target , double Position[3] )
{
  double    A[50] , Cp[50]  , sum[3] , T_break , T_seg , T_sub , Tc;
  int       i , j;
  long int  C , G , N , offset = 0;

  Time += TDBOFFSET; /****!!!!! Account for difference between TDB and UT**/
  /*--------------------------------------------------------------------------*/
  /* This function doesn't "do" nutations or librations.                      */
  /*--------------------------------------------------------------------------*/

  if ( Target >= 11 )             /* Also protects against weird input errors */
     {
       printf("\n This function does not compute nutations or librations.\n");
       return;
     }
 
  /*--------------------------------------------------------------------------*/
  /* Initialize local coefficient array.                                      */
  /*--------------------------------------------------------------------------*/

  for ( i=0 ; i<50 ; i++ )
      {
        A[i] = 0.0;
      }

  /*--------------------------------------------------------------------------*/
  /* Determine if a new record needs to be input (if so, get it).             */
  /*--------------------------------------------------------------------------*/
    
  if (Time < T_beg || Time > T_end)  Read_Coefficients(Time);

  /*--------------------------------------------------------------------------*/
  /* Read the coefficients from the binary record.                            */
  /*--------------------------------------------------------------------------*/
  
  C = R1.coeffPtr[Target][0] - 1;          /*   Coefficient array entry point */
  N = R1.coeffPtr[Target][1];              /* Number of coeff's per component */
  G = R1.coeffPtr[Target][2];              /*      Granules in current record */

  /*...................................................Debug print (optional) */

  if ( FALSE )
     {
       printf("\n  In: Interpolate_Position\n");
       printf("\n  Target = %2d",Target);
       printf("\n  C      = %4d (before)",C);
       printf("\n  N      = %4d",N);
       printf("\n  G      = %4d\n",G);
     }

  /*--------------------------------------------------------------------------*/
  /*  Compute the normalized time, then load the Tchebeyshev coefficients     */
  /*  into array A[]. If T_span is covered by a single granule this is easy.  */
  /*  If not, the granule that contains the interpolation time is found, and  */
  /*  an offset from the array entry point for the ephemeris body is used to  */
  /*  load the coefficients.                                                  */
  /*--------------------------------------------------------------------------*/

  if ( G == 1 )
     {
       Tc = 2.0*(Time - T_beg) / T_span - 1.0;
       for (i=C ; i<(C+3*N) ; i++)  A[i-C] = Coeff_Array[i];
     }
  else if ( G > 1 )
     {
       T_sub = T_span / ((double) G);          /* Compute subgranule interval */
       
       for ( j=G ; j>0 ; j-- ) 
           {
             T_break = T_beg + ((double) j-1) * T_sub;
             if ( Time > T_break ) 
                {
                  T_seg  = T_break;
                  offset = j-1;
                  break;
                }
            }
            
       Tc = 2.0*(Time - T_seg) / T_sub - 1.0;
       C  = C + 3 * offset * N;
       
       for (i=C ; i<(C+3*N) ; i++) A[i-C] = Coeff_Array[i];
     }
  else                                   /* Something has gone terribly wrong */
     {
       printf("\n Number of granules must be >= 1: check header data.\n");
     }

  /*...................................................Debug print (optional) */

  if ( FALSE )
     {
       printf("\n  C      = %4d (after)",C);
       printf("\n  offset = %4d",offset);
       printf("\n  Time   = %12.7f",Time);
       printf("\n  T_sub  = %12.7f",T_sub);
       printf("\n  T_seg  = %12.7f",T_seg);
       printf("\n  Tc     = %12.7f\n",Tc);
       printf("\n  Array Coefficients:\n");
       for ( i=0 ; i<3*N ; i++ )
           {
             printf("\n  A[%2d] = % 22.15e",i,A[i]);
           }
       printf("\n\n");
     }

  /*..........................................................................*/

  /*--------------------------------------------------------------------------*/
  /* Compute interpolated the position.                                       */
  /*--------------------------------------------------------------------------*/
  
  for ( i=0 ; i<3 ; i++ ) 
      {                           
        Cp[0]  = 1.0;                                 /* Begin polynomial sum */
        Cp[1]  = Tc;
        sum[i] = A[i*N] + A[1+i*N]*Tc;

        for ( j=2 ; j<N ; j++ )                                  /* Finish it */
            {
              Cp[j]  = 2.0 * Tc * Cp[j-1] - Cp[j-2];
              sum[i] = sum[i] + A[j+i*N] * Cp[j];
            }
        Position[i] = sum[i];
      }

  return;
}

void Interpolate_State(double Time, 
		       int Target, 
		       double Position[3], 
		       double Velocity[3])
{
  double    A[50]   , B[50] , Cp[50] , P_Sum[3] , V_Sum[3] , Up[50] ,
            T_break , T_seg , T_sub  , Tc;
  int       i , j;
  long int  C , G , N , offset = 0;

  Time += TDBOFFSET; /****!!!!! Account for difference between TDB and UT**/

  /*--------------------------------------------------------------------------*/
  /* This function doesn't "do" nutations or librations.                      */
  /*--------------------------------------------------------------------------*/

  if ( Target >= 11 )             /* Also protects against weird input errors */
     {
       printf("\n This function does not compute nutations or librations.\n");
       return;
     }

  /*--------------------------------------------------------------------------*/
  /* Initialize local coefficient array.                                      */
  /*--------------------------------------------------------------------------*/

  for ( i=0 ; i<50 ; i++ )
      {
        A[i] = 0.0;
        B[i] = 0.0;
      }

  /*--------------------------------------------------------------------------*/
  /* Determine if a new record needs to be input.                             */
  /*--------------------------------------------------------------------------*/
  
  if (Time < T_beg || Time > T_end)  Read_Coefficients(Time);

  /*--------------------------------------------------------------------------*/
  /* Read the coefficients from the binary record.                            */
  /*--------------------------------------------------------------------------*/
  
  C = R1.coeffPtr[Target][0] - 1;               /*    Coeff array entry point */
  N = R1.coeffPtr[Target][1];                   /*          Number of coeff's */
  G = R1.coeffPtr[Target][2];                   /* Granules in current record */

  /*...................................................Debug print (optional) */

  if ( FALSE )
     {
       printf("\n  In: Interpolate_State\n");
       printf("\n  Target = %2d",Target);
       printf("\n  C      = %4d (before)",C);
       printf("\n  N      = %4d",N);
       printf("\n  G      = %4d\n",G);
     }

  /*--------------------------------------------------------------------------*/
  /*  Compute the normalized time, then load the Tchebeyshev coefficients     */
  /*  into array A[]. If T_span is covered by a single granule this is easy.  */
  /*  If not, the granule that contains the interpolation time is found, and  */
  /*  an offset from the array entry point for the ephemeris body is used to  */
  /*  load the coefficients.                                                  */
  /*--------------------------------------------------------------------------*/

  if ( G == 1 )
     {
       Tc = 2.0*(Time - T_beg) / T_span - 1.0;
       for (i=C ; i<(C+3*N) ; i++)  A[i-C] = Coeff_Array[i];
     }
  else if ( G > 1 )
     {
       T_sub = T_span / ((double) G);          /* Compute subgranule interval */
       
       for ( j=G ; j>0 ; j-- ) 
           {
             T_break = T_beg + ((double) j-1) * T_sub;
             if ( Time > T_break ) 
                {
                  T_seg  = T_break;
                  offset = j-1;
                  break;
                }
            }
            
       Tc = 2.0*(Time - T_seg) / T_sub - 1.0;
       C  = C + 3 * offset * N;
       
       for (i=C ; i<(C+3*N) ; i++) A[i-C] = Coeff_Array[i];
     }
  else                                   /* Something has gone terribly wrong */
     {
       printf("\n Number of granules must be >= 1: check header data.\n");
     }

  /*...................................................Debug print (optional) */

  if ( FALSE )
     {
       printf("\n  C      = %4d (after)",C);
       printf("\n  offset = %4d",offset);
       printf("\n  Time   = %12.7f",Time);
       printf("\n  T_sub  = %12.7f",T_sub);
       printf("\n  T_seg  = %12.7f",T_seg);
       printf("\n  Tc     = %12.7f\n",Tc);
       printf("\n  Array Coefficients:\n");
       for ( i=0 ; i<3*N ; i++ )
           {
             printf("\n  A[%2d] = % 22.15e",i,A[i]);
           }
       printf("\n\n");
     }

  /*..........................................................................*/

  /*--------------------------------------------------------------------------*/
  /* Compute the interpolated position & velocity                             */
  /*--------------------------------------------------------------------------*/
  
  for ( i=0 ; i<3 ; i++ )                /* Compute interpolating polynomials */
      {
        Cp[0] = 1.0;           
        Cp[1] = Tc;
        Cp[2] = 2.0 * Tc*Tc - 1.0;
        
        Up[0] = 0.0;
        Up[1] = 1.0;
        Up[2] = 4.0 * Tc;

        for ( j=3 ; j<N ; j++ )
            {
              Cp[j] = 2.0 * Tc * Cp[j-1] - Cp[j-2];
              Up[j] = 2.0 * Tc * Up[j-1] + 2.0 * Cp[j-1] - Up[j-2];
            }

        P_Sum[i] = 0.0;           /* Compute interpolated position & velocity */
        V_Sum[i] = 0.0;

        for ( j=N-1 ; j>-1 ; j-- )  P_Sum[i] = P_Sum[i] + A[j+i*N] * Cp[j];
        for ( j=N-1 ; j>0  ; j-- )  V_Sum[i] = V_Sum[i] + A[j+i*N] * Up[j];

        Position[i] = P_Sum[i];
        Velocity[i] = V_Sum[i] * 2.0 * ((double) G) / (T_span * 86400.0);
      }

  /*--------------------------------------------------------------------------*/
  /*  Return computed values.                                                 */
  /*--------------------------------------------------------------------------*/

  return;
}

/********************************************************* END: ephem_read.c **/

/* Give the Earth geocenter wrt SSBary in AU. */
void
geocenter_ssbary(double jd,
		double *xyz)
{
  double embary[3],moon[3];
  static int init=0;
  int i;

  if (!init) {
    if (Initialize_Ephemeris()) exit(1);
    init = 1;
  }

  Interpolate_Position(jd, EARTH, embary);
  Interpolate_Position(jd, MOON, moon);
  for (i=0; i<3; i++) {
    xyz[i] = (embary[i] - moon[i]/(1.+R1.EMRAT)) / R1.AU;
  }

  return;
}

/********** Stealing some routines from skycalc to do the shift
******** from geocenter to observatory.
********/
#define  J2000             2451545.        /* Julian date at standard epoch */
#define  SEC_IN_DAY        86400.
#define  FLATTEN           0.003352813   /* flattening of earth, 1/298.257 */
#define  EQUAT_RAD         6378137.    /* equatorial radius of earth, meters */

double 
lst(double jd,
    double longit)
{
	/* returns the local MEAN sidereal time (dec hrs) at julian date jd
	   at west longitude long (decimal hours).  Follows
	   definitions in 1992 Astronomical Almanac, pp. B7 and L2.
	   Expression for GMST at 0h ut referenced to Aoki et al, A&A 105,
	   p.359, 1982.  On workstations, accuracy (numerical only!)
	   is about a millisecond in the 1990s. */

	double t, ut, jdmid, jdint, jdfrac, sid_g, sid;
	long jdin, sid_int;

	jdin = jd;         /* fossil code from earlier package which
			split jd into integer and fractional parts ... */
	jdint = jdin;
	jdfrac = jd - jdint;
	if(jdfrac < 0.5) {
		jdmid = jdint - 0.5;
		ut = jdfrac + 0.5;
	}
	else {
		jdmid = jdint + 0.5;
		ut = jdfrac - 0.5;
	}
	t = (jdmid - J2000)/36525;
	sid_g = (24110.54841+8640184.812866*t+0.093104*t*t-6.2e-6*t*t*t)/SEC_IN_DAY;
	sid_int = sid_g;
	sid_g = sid_g - (double) sid_int;
	sid_g = sid_g + 1.0027379093 * ut - longit/24.;
	sid_int = sid_g;
	sid_g = (sid_g - (double) sid_int) * 24.;
	if(sid_g < 0.) sid_g = sid_g + 24.;
	return(sid_g);
}

void 
topo(double geolong, double geolat, double height,
	double *x_geo, double *y_geo, double *z_geo)
/* computes the geocentric coordinates from the geodetic
(standard map-type) longitude, latitude, and height.
These are assumed to be in decimal hours, decimal degrees, and
meters respectively.  Notation generally follows 1992 Astr Almanac,
p. K11 */

/**** I've changed this to assume long & lat in radians, and
***** to return distances in AU.  
***** And note that this gives proper equatorial vector only if lon
***** is really LMST measured positive eastward.
****/
{

	double denom, C_geo, S_geo;

	denom = (1. - FLATTEN) * sin(geolat);
	denom = cos(geolat) * cos(geolat) + denom*denom;
	C_geo = 1. / sqrt(denom);
	S_geo = (1. - FLATTEN) * (1. - FLATTEN) * C_geo;
	C_geo = C_geo* EQUAT_RAD + height ;  /* deviation from almanac
		       notation -- include height here. */
	S_geo = S_geo* EQUAT_RAD + height ;
	*x_geo = C_geo * cos(geolat) * cos(geolong);
	*y_geo = C_geo * cos(geolat) * sin(geolong);
	*z_geo = S_geo * sin(geolat);

	/* convert to AU, keeping in mind that Horizons was km and this is m*/
	*x_geo /= (1000.*R1.AU);
	*y_geo /= (1000.*R1.AU);
	*z_geo /= (1000.*R1.AU);
}

/* OK, now put all of the above together to give the observatory coordinates
 * wrt SSBary in ICRS/AU.
 */
void
earth_ssbary(double jd,
	     int    obscode,
	     double *xtopo, double *ytopo, double *ztopo)
{
  double	xgeo[3], xobs,yobs,zobs;

  /* Return zero if we are "observing" from ssbary */
  if (obscode==OBSCODE_SSBARY) {
    *xtopo = *ytopo = *ztopo = 0.;
    return;
  }

  /* Use Horizons to get the coordinates of geocenter */
  geocenter_ssbary(jd,xgeo);

  /*fprintf(stderr,"geocenter ICRS:  %g %g %g\n",xgeo[0],xgeo[1],xgeo[2]);*/
  /* Now get the geocenter->observatory vector */
  observatory_geocenter(jd, obscode, &xobs, &yobs, &zobs);

  /*fprintf(stderr,"observ ICRS:  %g %g %g\n", xobs, yobs, zobs);*/

  /* Add in to give observatory coords in ICRS */
  *xtopo = xgeo[0] + xobs;
  *ytopo = xgeo[1] + yobs;
  *ztopo = xgeo[2] + zobs;

  return;
}

/* return vector from geocenter->observatory in ICRS coordinates*/
void
observatory_geocenter(double jd,
		      int obscode,
		      double *xobs,
		      double *yobs,
		      double *zobs) {

  if (nsites<=0 && nspacecraft<=0) read_observatories(NULL); 

  if (obscode < OBSCODE_ORBITAL) {
    /* This is a ground-based observatory */
    double	obslat, obslon, obsalt, obslmst;	/*observatory info*/
    int i;

    if (obscode==OBSCODE_GEOCENTER) {
      *xobs=*yobs=*zobs=0.;
      return;
    }

    /* get the lat, lon, & altitude of this observatory */
    for (i=0; obscode!=sitelist[i].code && i<nsites; i++)  ;
    if (i>=nsites) {
      fprintf(stderr,"Unknown observatory code %d\n",obscode);
      exit(1);
    }
    obslon = sitelist[i].lon;
    obslat = sitelist[i].lat;
    obsalt = sitelist[i].altitude;

    /* Get the LMST and calculate to ICRS vector */
    obslmst=lst(jd,obslon);
    topo(obslmst*PI/12., obslat, obsalt, xobs, yobs, zobs);

  } else {
    /* Orbiting observatory */
    double pole, rapole, phase;
    int i;
    SPACECRAFT *s;
    for (i=0; obscode!=spacecraftlist[i].code && i<nspacecraft; i++)  ;
    if (i>=nspacecraft) {
      fprintf(stderr,"Unknown spacecraft code %d\n",obscode);
      exit(1);
    }
    s = &(spacecraftlist[i]);
    pole = s->i;
    phase = -(jd - s->jd0) / s->precess;/* RA of pole DEcreases w/time 
					 * for precess > 0 */
    phase -= floor(phase);
    rapole = s->ra0 + phase*TPI;

    phase = (jd - s->jd0) / s->P;	/* orbit is increasing RA 
					 * for P > 0 */
    phase -= floor(phase);
    phase *= TPI;

    *xobs = s->a * (cos(phase)*cos(rapole) 
		    - cos(pole)*sin(phase)*sin(rapole));
    *yobs = s->a * (cos(phase)*sin(rapole) 
		    + cos(pole)*sin(phase)*cos(rapole));
    *zobs = s->a * sin(phase)*sin(pole);
    /*fprintf(stderr,"Obs. Posn at jd %.5f : phase %.1f rapole %.1f pole %.1f \n",
      jd, phase/DTOR, rapole/DTOR, pole/DTOR);*/
  /*fprintf(stderr," x,y,z: %g %g %g\n", *xobs, *yobs, *zobs);*/
  }

  return;
}


/* Read the look-up table for observatories.  
 * Note that ground-based longitudes are stored in hours, 
 * altitude in meters, lat in radians */
void
read_observatories(char *fname)
{
  int i;
  FILE *sitefile;
  char  inbuff[BUFFSIZE], latstring[40], lonstring[40];
  int nchar, obscode;
  char fileName[FNAMESIZE];

  extern double dmsdeg(char *string);
  extern double hmsdeg(char *string);

  nsites = nspacecraft = 0;

  /** use passed filename, or a previously specified filename,
   ** or environment-specified file, or the default filename, in that order
   */
  if (fname != NULL)
    strncpy(fileName, fname, FNAMESIZE-1);
  else if (strlen(observatory_file)>0)
    strncpy(fileName, observatory_file, FNAMESIZE-1);
  else if (getenv(OBS_ENVIRON)!=NULL)
    strncpy(fileName, getenv(OBS_ENVIRON), FNAMESIZE-1);
  else
    strncpy(fileName, DEFAULT_OBSERVATORY_FILE, FNAMESIZE-1);
  
  fileName[FNAMESIZE-1]=0;

  if ((sitefile = fopen(fileName,"r"))==NULL) {
    fprintf(stderr,"Error opening observatories file %s\n",fileName);
    exit(1);
  }

  while (fgets_nocomment(inbuff, BUFFSIZE-1, sitefile, NULL)!=NULL) {
    if (sscanf(inbuff, "%d", &obscode)!=1) {
      fprintf(stderr,"Bad line in obseratories file %s:\n->%s\n",
	      fileName, inbuff);
      exit(1);
    }
    if (obscode < OBSCODE_ORBITAL) {
      /* This line specifies a ground-based location */
      SITE *sss;
      sss = &(sitelist[nsites]);
      if (sscanf(inbuff,"%d %s %s %lf %n", &(sss->code), lonstring, latstring,
		 &(sss->altitude), &nchar)!=4) {
	fprintf(stderr,"Bad line in obseratories file %s:\n->%s\n",
		fileName, inbuff);
	exit(1);
      }
      if (nchar<strlen(inbuff)) {
	int nn;
	nn = strcspn(inbuff+nchar,"\n");
	if (nn>=80) nn=79;
	strncpy(sss->name, inbuff+nchar, nn);
	sss->name[nn]=0;
      } else {
	sss->name[0] = 0;
      }
      sss->lon = dmsdeg(lonstring)/15.;	/*convert longitude to hours*/
      sss->lat = dmsdeg(latstring)*DTOR;
      nsites++;

    } else {
      /* This is an orbiting observatory */
      SPACECRAFT *sss;
      sss = &(spacecraftlist[nspacecraft]);
      if (sscanf(inbuff,"%d %s %lf %lf %lf %s %n", 
		 &(sss->code), lonstring, &(sss->P),
		 &(sss->precess), &(sss->jd0), latstring, 
		 &nchar)!=6) {
	fprintf(stderr,"Bad line in obseratories file %s:\n->%s\n",
		fileName, inbuff);
	exit(1);
      }
      if (nchar<strlen(inbuff)) {
	int nn;
	nn = strcspn(inbuff+nchar,"\n");
	if (nn>=80) nn=79;
	strncpy(sss->name, inbuff+nchar, nn);
	sss->name[nn]=0;
      } else {
	sss->name[0] = 0;
      }
      sss->i = dmsdeg(lonstring)*DTOR;	/*convert to radians*/
      sss->ra0 = dmsdeg(latstring)*DTOR;

      /* calculate the semi-major axis, in AU */
      sss->a = pow( pow(sss->P*DAY,2.) * GM * EARTHMASS / (TPI*TPI),
		    1./3.);
      nspacecraft++;
    }
  }
  fclose(sitefile);
  
  if (nsites<1 && nspacecraft<1) {
    fprintf(stderr,"Error: no observatory sites found\n");
    exit(1);
  }
}

/* Get arbitrary planetary barycenter position.  Get the
velocity (in AU/YR) too, if desired. */
void
bodycenter_ssbary(double jd,
		  double *xyz,
		  int body,
		  double *vxyz)
{
  double posn[3], vel[3];
  static int init=0;
  int i;

  if (!init) {
    if (Initialize_Ephemeris()) exit(1);
    init = 1;
  }

  if (vxyz==NULL) {
    Interpolate_Position(jd, body, posn);
  } else {
    Interpolate_State(jd, body, posn, vel);
  }
  for (i=0; i<3; i++) {
    xyz[i] = posn[i] / R1.AU;
  }
  if (vxyz!=NULL) /* convert km/s to AU/YR: */
    for (i=0; i<3; i++) {
      vxyz[i] = vel[i] * (86400. * 365.25) / R1.AU ;
    }
  return;
}

/* Return the angle from zenith to the horizon for this observatory*/
double
zenith_horizon(int obscode) {
  if (obscode < OBSCODE_ORBITAL)
    return PI/2.;
  else {
    /* Do spherical-Earth limb calculation for orbiting spacecraft */
    double a;
    int i;
    if (nsites<=0 && nspacecraft<=0) read_observatories(NULL); 
    for (i=0; obscode!=spacecraftlist[i].code && i<nspacecraft; i++)  ;
    if (i>=nspacecraft) {
      fprintf(stderr,"Unknown spacecraft code %d\n",obscode);
      exit(1);
    }
    a = spacecraftlist[i].a;
    a /= EQUAT_RAD/(1000.*R1.AU);
    if (a < 1.) {
      fprintf(stderr,"Your spacecraft is underground.  Oops.\n");
      exit(1);
    }
    return PI - asin(1./a);
  }
}
