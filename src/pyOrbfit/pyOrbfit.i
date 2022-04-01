/* File : pyOrbfit.i */
%module pyOrbfit
%include typemaps.i
%{
#define SWIG_FILE_WITH_INIT
#include "orbfit.h"
%}
%include "numpy.i"
%include <carrays.i>
%init %{
import_array();
%}
%array_class(OBSERVATION, OBSERVATION_ARRAY);
%array_class(double, doubleArray);
%define %apply_numpy_typemaps(TYPE)

%apply (TYPE IN_ARRAY2[ANY][ANY]) {(TYPE matrix[ANY][ANY])};
%apply (TYPE* IN_ARRAY2, int DIM1, int DIM2) {(TYPE* matrix, int rows, int cols)};
%apply (int DIM1, int DIM2, TYPE* IN_ARRAY2) {(int rows, int cols, TYPE* matrix)};

%apply (TYPE INPLACE_ARRAY2[ANY][ANY]) {(TYPE array[6][6])};
%apply (TYPE* INPLACE_ARRAY2, int DIM1, int DIM2) {(TYPE* array, int rows, int cols)};
%apply (int DIM1, int DIM2, TYPE* INPLACE_ARRAY2) {(int rows, int cols, TYPE* array)};

%enddef    /* %apply_numpy_typemaps() macro */
%apply_numpy_typemaps(signed char       )
%apply_numpy_typemaps(unsigned char     )
%apply_numpy_typemaps(short             )
%apply_numpy_typemaps(unsigned short    )
%apply_numpy_typemaps(int               )
%apply_numpy_typemaps(unsigned int      )
%apply_numpy_typemaps(long              )
%apply_numpy_typemaps(unsigned long     )
%apply_numpy_typemaps(long long         )
%apply_numpy_typemaps(unsigned long long)
%apply_numpy_typemaps(float             )
%apply_numpy_typemaps(double            )
%apply int *OUTPUT {int *nobservations, int *ndof, int *dof}; 
%apply double *OUTPUT {double *chisqfit, double *chisq};
%apply double *OUTPUT {double *x, double *y};
%apply double *OUTPUT {double *lat_ec,  double *lon_ec};
%apply double *OUTPUT {double *ra_eq,  double *dec_eq};
%apply int *OUTPUT {int *nobs};
%apply char *OUTPUT {char *outbuff};
/* Put headers and other declarations here */
%{
  extern double xBary, yBary, zBary;
  extern double lat0, lon0, jd0; 
  extern void set_observatory_file(char *fname);
  extern void set_ephem_file(char *fname);
  extern double **dmatrix(int nrl,int nrh,int ncl, int nch);
  extern double *dvector(int nl,int nh);
  extern void fit_radec(char *fname, int *nobservations, double *chisqfit, int *ndof, PBASIS *pbasis, ORBIT *orbit);
  extern int read_radec(OBSERVATION_ARRAY obsarray[], char *fname, int *nobs);
  extern int fit_observations(OBSERVATION_ARRAY obsarray[], int nobs, PBASIS *p, double **covar, double *chisq, int *dof, FILE *logfile);
  extern void pbasis_to_bary(PBASIS *p, XVBASIS *xv, double **partials);
  extern void flatten_cov(double **cov, int ndim, double *cov1d);
  extern void unflatten_cov(double *cov1d, int ndim, double **cov);
  extern void add_to_obsarray(OBSERVATION_ARRAY obsarray[], int iobs, OBSERVATION obs);
  extern void kbo2d(PBASIS *pin, OBSERVATION *obs, double *x, double dx[], double *y, double dy[]);
  extern void deghms(double degr, char *outbuff);
  extern void orbitElements(XVBASIS *xv, ORBIT  *orb);
%}

/* Parse the header file to generate wrappers */
/*%apply (double* INPLACE_ARRAY1, int DIM1) {(double* cov1d, int ndim)}; */
%include "orbfit.h";
extern double xBary, yBary, zBary;
extern double lat0, lon0, jd0; 
extern void set_ephem_file(char *fname);
extern void set_observatory_file(char *fname);
extern double **dmatrix(int nrl,int nrh,int ncl, int nch);
extern double *dvector(int nl, int nh);
extern void fit_radec(char *fname, int *nobservations, double *chisqfit, int *ndof, PBASIS *pbasis, ORBIT *orbit);
extern int read_radec(OBSERVATION_ARRAY obsarray[], char *fname, int *nobs);
extern int fit_observations(OBSERVATION_ARRAY obsarray[], int nobs, PBASIS *p, double **covar, double *chisq, int *dof, FILE *logfile);
extern void pbasis_to_bary(PBASIS *p, XVBASIS *xv, double **partials);
extern void flatten_cov(double **cov, int ndim, double *cov1d);
extern void unflatten_cov(double *cov1d, int ndim, double **cov);
extern void add_to_obsarray(OBSERVATION_ARRAY obsarray[], int iobs, OBSERVATION obs);
extern void kbo2d(PBASIS *pin, OBSERVATION *obs, double *x, double dx[], double *y, double dy[]);
extern void deghms(double degr, char *outbuff);
extern void orbitElements(XVBASIS *xv,ORBIT  *orb);
