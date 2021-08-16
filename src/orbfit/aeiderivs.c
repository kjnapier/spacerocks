/* 	$Id: aeiderivs.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: aeiderivs.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* Subroutine to determine partial derivative matrix of orbital params
 * w.r.t. barycentric phase space params.  Incredible algebra by Bharat.
 */

#include "orbfit.h"

void
aei_derivs( XVBASIS *xv,
	    double **daei_dxv)
{
  double x,y,z,xdot,ydot,zdot;
  double mu=GM*SSMASS; /* AU^3/yr^2 */
  /* adjusted by factor in parentheses to scale to entire SS mass from Sun*/

  x = xv->x;
  y = xv->y;
  z = xv->z;
  xdot = xv->xdot;
  ydot = xv->ydot;
  zdot = xv->zdot;
  /* x,y,z,xdot,ydot,zdot are heliocentric eq. x,y,z in AU 
     and xdot,ydot,zdot in AU/yr */

  /* ax is del a by del x */
  daei_dxv[1][1] = 
    (2*x)/(pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*
	   pow(2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
	       (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,2));
       
  daei_dxv[1][2] = 
  (2*y)/(pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*
	    pow(2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
		(pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,2));
       
  daei_dxv[1][3] = 
  (2*z)/(pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*
	    pow(2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
		(pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,2));
       
  daei_dxv[1][4] = 
  (2*xdot)/(mu*pow(2/
			 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
			 (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,2));
       
  daei_dxv[1][5] = 
  (2*ydot)/(mu*pow(2/
			 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
			 (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,2));
       
  daei_dxv[1][6] = 
    (2*zdot)/(mu*pow(2/
		     sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
		     (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,2));
       
       
/* Partials for e now */

  daei_dxv[2][1] = 
    (2*(-(pow(xdot,2)/mu) + 
        pow(x,2)/
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - 
        1/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + 
        (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)*
      (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
        x*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
       + 2*(-((xdot*ydot)/mu) + 
        (x*y)/pow(pow(x,2) + pow(y,2) + pow(z,2),1.5))*
      (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
        y*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
       + 2*((x*z)/
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - 
        (xdot*zdot)/mu)*(-((zdot*(x*xdot + y*ydot + z*zdot))/
           mu) + z*(-(1/
              sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
     )/(2.*sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
            mu) + x*(-(1/
               sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2) + pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
         y*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2) + pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
         z*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2)));
	 	 
  daei_dxv[2][2] = 
    (2*(-((xdot*ydot)/mu) + (x*y)/
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5))*
      (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
        x*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
       + 2*(-(pow(ydot,2)/mu) + 
        pow(y,2)/
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - 
        1/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + 
        (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)*
      (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
        y*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
       + 2*((y*z)/
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - 
        (ydot*zdot)/mu)*(-((zdot*(x*xdot + y*ydot + z*zdot))/
           mu) + z*(-(1/
              sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
     )/(2.*sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
            mu) + x*(-(1/
               sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2) + pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
         y*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2) + pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
         z*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2)));
	 
  daei_dxv[2][3] = 
    (2*((x*z)/pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - 
        (xdot*zdot)/mu)*(-((xdot*(x*xdot + y*ydot + z*zdot))/
           mu) + x*(-(1/
              sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
       + 2*((y*z)/
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - 
        (ydot*zdot)/mu)*(-((ydot*(x*xdot + y*ydot + z*zdot))/
           mu) + y*(-(1/
              sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
       + 2*(pow(z,2)/
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - 
        1/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
        pow(zdot,2)/mu + 
        (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)*
      (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
        z*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
     )/(2.*sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
            mu) + x*(-(1/
               sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2) + pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
         y*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2) + pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
         z*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2)));
	 
  daei_dxv[2][4] = 
    (2*((x*xdot)/mu - (x*xdot + y*ydot + z*zdot)/mu)*
      (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
        x*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
       + 2*((2*xdot*y)/mu - (x*ydot)/mu)*
      (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
        y*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
       + 2*((2*xdot*z)/mu - (x*zdot)/mu)*
      (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
        z*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
     )/(2.*sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
            mu) + x*(-(1/
               sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2) + pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
         y*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2) + pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
         z*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2)));
	 
  daei_dxv[2][5] = 
    (2*(-((xdot*y)/mu) + (2*x*ydot)/mu)*
      (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
        x*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
       + 2*((y*ydot)/mu - (x*xdot + y*ydot + z*zdot)/mu)*
      (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
        y*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
       + 2*((2*ydot*z)/mu - (y*zdot)/mu)*
      (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
        z*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
     )/(2.*sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
            mu) + x*(-(1/
               sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2) + pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
         y*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2) + pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
         z*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2)));
	 
  daei_dxv[2][6] = 
    (2*(-((xdot*z)/mu) + (2*x*zdot)/mu)*
      (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
        x*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
       + 2*(-((ydot*z)/mu) + (2*y*zdot)/mu)*
      (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
        y*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
       + 2*((z*zdot)/mu - (x*xdot + y*ydot + z*zdot)/mu)*
      (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
        z*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
           (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu))
     )/(2.*sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
            mu) + x*(-(1/
               sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2) + pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
         y*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2) + pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
         z*(-(1/sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + 
            (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu)
         ,2)));
	 
/* Partials of i now */

  daei_dxv[3][1] = 
    -((-((-(xdot*y) + x*ydot)*
           (2*ydot*(-(xdot*y) + x*ydot) - 
             2*zdot*(xdot*z - x*zdot)))/
        (2.*pow(pow(-(xdot*y) + x*ydot,2) + 
            pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2),1.5)) + 
       ydot/sqrt(pow(-(xdot*y) + x*ydot,2) + 
          pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
          ))/
     sqrt(1 - pow(-(xdot*y) + x*ydot,2)/
        (pow(-(xdot*y) + x*ydot,2) + 
          pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
          )));
	  
  daei_dxv[3][2] = 
    -((-((-(xdot*y) + x*ydot)*
           (-2*xdot*(-(xdot*y) + x*ydot) + 
             2*zdot*(-(ydot*z) + y*zdot)))/
        (2.*pow(pow(-(xdot*y) + x*ydot,2) + 
            pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2),1.5)) - 
       xdot/sqrt(pow(-(xdot*y) + x*ydot,2) + 
          pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
          ))/
     sqrt(1 - pow(-(xdot*y) + x*ydot,2)/
        (pow(-(xdot*y) + x*ydot,2) + 
          pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
          )));
	  
  daei_dxv[3][3] = 
    ((-(xdot*y) + x*ydot)*(2*xdot*(xdot*z - x*zdot) - 
       2*ydot*(-(ydot*z) + y*zdot)))/
   (2.*pow(pow(-(xdot*y) + x*ydot,2) + 
       pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2),
      1.5)*sqrt(1 - pow(-(xdot*y) + x*ydot,2)/
        (pow(-(xdot*y) + x*ydot,2) + 
          pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
          )));
	  
  daei_dxv[3][4] = 
    -((-((-(xdot*y) + x*ydot)*
           (-2*y*(-(xdot*y) + x*ydot) + 2*z*(xdot*z - x*zdot)))/
        (2.*pow(pow(-(xdot*y) + x*ydot,2) + 
            pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2),1.5)) - 
       y/sqrt(pow(-(xdot*y) + x*ydot,2) + 
          pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
          ))/
     sqrt(1 - pow(-(xdot*y) + x*ydot,2)/
        (pow(-(xdot*y) + x*ydot,2) + 
          pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
          )));
	  
  daei_dxv[3][5] = 
    -((-((-(xdot*y) + x*ydot)*
           (2*x*(-(xdot*y) + x*ydot) - 2*z*(-(ydot*z) + y*zdot))
           )/
        (2.*pow(pow(-(xdot*y) + x*ydot,2) + 
            pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2),1.5)) + 
       x/sqrt(pow(-(xdot*y) + x*ydot,2) + 
          pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
          ))/
     sqrt(1 - pow(-(xdot*y) + x*ydot,2)/
        (pow(-(xdot*y) + x*ydot,2) + 
          pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
          )));
	  
  daei_dxv[3][6] = 
    ((-(xdot*y) + x*ydot)*(-2*x*(xdot*z - x*zdot) + 
       2*y*(-(ydot*z) + y*zdot)))/
   (2.*pow(pow(-(xdot*y) + x*ydot,2) + 
       pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2),
      1.5)*sqrt(1 - pow(-(xdot*y) + x*ydot,2)/
        (pow(-(xdot*y) + x*ydot,2) + 
          pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
          )));
	  

/* Partials of capital Omega (long of asc node) now */

  daei_dxv[4][1] = 
    -(((zdot*(xdot*z - x*zdot)*(-(xdot*z) + x*zdot))/
        pow(pow(xdot*z - x*zdot,2) + 
          pow(-(ydot*z) + y*zdot,2),1.5) + 
       zdot/sqrt(pow(xdot*z - x*zdot,2) + 
          pow(-(ydot*z) + y*zdot,2)))/
     sqrt(1 - pow(-(xdot*z) + x*zdot,2)/
        (pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2))
       ));
       
  daei_dxv[4][2] = 
    (zdot*(-(xdot*z) + x*zdot)*(-(ydot*z) + y*zdot))/
   (pow(pow(xdot*z - x*zdot,2) + 
       pow(-(ydot*z) + y*zdot,2),1.5)*
     sqrt(1 - pow(-(xdot*z) + x*zdot,2)/
        (pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2))
       ));
       
  daei_dxv[4][3] = 
    -((-((-(xdot*z) + x*zdot)*
           (2*xdot*(xdot*z - x*zdot) - 
             2*ydot*(-(ydot*z) + y*zdot)))/
        (2.*pow(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2),1.5)) - 
       xdot/sqrt(pow(xdot*z - x*zdot,2) + 
          pow(-(ydot*z) + y*zdot,2)))/
     sqrt(1 - pow(-(xdot*z) + x*zdot,2)/
        (pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2))
       ));
       
  daei_dxv[4][4] = 
    -((-((z*(xdot*z - x*zdot)*(-(xdot*z) + x*zdot))/
          pow(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2),1.5)) - 
       z/sqrt(pow(xdot*z - x*zdot,2) + 
          pow(-(ydot*z) + y*zdot,2)))/
     sqrt(1 - pow(-(xdot*z) + x*zdot,2)/
        (pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2))
       ));
       
  daei_dxv[4][5] = 
    -((z*(-(xdot*z) + x*zdot)*(-(ydot*z) + y*zdot))/
     (pow(pow(xdot*z - x*zdot,2) + 
         pow(-(ydot*z) + y*zdot,2),1.5)*
       sqrt(1 - pow(-(xdot*z) + x*zdot,2)/
          (pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2)))));
	    
  daei_dxv[4][6] = 
    -((-((-(xdot*z) + x*zdot)*
           (-2*x*(xdot*z - x*zdot) + 2*y*(-(ydot*z) + y*zdot)))/
        (2.*pow(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2),1.5)) + 
       x/sqrt(pow(xdot*z - x*zdot,2) + 
          pow(-(ydot*z) + y*zdot,2)))/
     sqrt(1 - pow(-(xdot*z) + x*zdot,2)/
        (pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2))
       ));
       

/* partials of small omega (arg of per) now */

  daei_dxv[5][1] = 
    -((-(((-(xdot*z) + x*zdot)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             (-(ydot*z) + y*zdot)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)))*
           (2*(-(pow(xdot,2)/mu) + 
                pow(x,2)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - 1/
                 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             2*(-((xdot*ydot)/mu) + 
                (x*y)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5))*(-((ydot*(x*xdot + y*ydot + z*zdot))/
                   mu) + y*
                 (-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             2*((x*z)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - (xdot*zdot)/mu)*
              (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu))))/
        (2.*sqrt(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          pow(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                 mu) + x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2),1.5)) + 
       ((-((xdot*ydot)/mu) + 
             (x*y)/
              pow(pow(x,2) + pow(y,2) + pow(z,2),1.5))*
           (-(ydot*z) + y*zdot) + 
          (-(xdot*z) + x*zdot)*
           (-(pow(xdot,2)/mu) + 
             pow(x,2)/
              pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)\
              - 1/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ) + zdot*(-((xdot*(x*xdot + y*ydot + z*zdot))/
                mu) + x*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)))/
        (sqrt(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2))) + 
       (zdot*(xdot*z - x*zdot)*
          ((-(xdot*z) + x*zdot)*
             (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu)) + 
            (-(ydot*z) + y*zdot)*
             (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu))))/
        (pow(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2),1.5)*
          sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2))))/
     sqrt(1 - pow((-(xdot*z) + x*zdot)*
           (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
             x*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)) + 
          (-(ydot*z) + y*zdot)*
           (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
             y*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)),2)/
        ((pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          (pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2)))));
	      
  daei_dxv[5][2] = 
    -((-(((-(xdot*z) + x*zdot)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             (-(ydot*z) + y*zdot)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)))*
           (2*(-((xdot*ydot)/mu) + 
                (x*y)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5))*(-((xdot*(x*xdot + y*ydot + z*zdot))/
                   mu) + x*
                 (-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             2*(-(pow(ydot,2)/mu) + 
                pow(y,2)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - 1/
                 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             2*((y*z)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - (ydot*zdot)/mu)*
              (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu))))/
        (2.*sqrt(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          pow(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                 mu) + x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2),1.5)) + 
       ((-((xdot*ydot)/mu) + 
             (x*y)/
              pow(pow(x,2) + pow(y,2) + pow(z,2),1.5))*
           (-(xdot*z) + x*zdot) + 
          (-(ydot*z) + y*zdot)*
           (-(pow(ydot,2)/mu) + 
             pow(y,2)/
              pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)\
              - 1/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ) + zdot*(-((ydot*(x*xdot + y*ydot + z*zdot))/
                mu) + y*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)))/
        (sqrt(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2))) - 
       (zdot*(-(ydot*z) + y*zdot)*
          ((-(xdot*z) + x*zdot)*
             (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu)) + 
            (-(ydot*z) + y*zdot)*
             (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu))))/
        (pow(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2),1.5)*
          sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2))))/
     sqrt(1 - pow((-(xdot*z) + x*zdot)*
           (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
             x*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)) + 
          (-(ydot*z) + y*zdot)*
           (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
             y*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)),2)/
        ((pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          (pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2)))));
	      
  daei_dxv[5][3] = 
    -((-(((-(xdot*z) + x*zdot)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             (-(ydot*z) + y*zdot)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)))*
           (2*((x*z)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - (xdot*zdot)/mu)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             2*((y*z)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - (ydot*zdot)/mu)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             2*(pow(z,2)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - 1/
                 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                pow(zdot,2)/mu + 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)*
              (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu))))/
        (2.*sqrt(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          pow(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                 mu) + x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2),1.5)) + 
       ((-(xdot*z) + x*zdot)*
           ((x*z)/
              pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)\
              - (xdot*zdot)/mu) + 
          (-(ydot*z) + y*zdot)*
           ((y*z)/
              pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)\
              - (ydot*zdot)/mu) - 
          xdot*(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
             x*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)) - 
          ydot*(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
             y*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)))/
        (sqrt(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2))) - 
       ((2*xdot*(xdot*z - x*zdot) - 
            2*ydot*(-(ydot*z) + y*zdot))*
          ((-(xdot*z) + x*zdot)*
             (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu)) + 
            (-(ydot*z) + y*zdot)*
             (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu))))/
        (2.*pow(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2),1.5)*
          sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2))))/
     sqrt(1 - pow((-(xdot*z) + x*zdot)*
           (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
             x*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)) + 
          (-(ydot*z) + y*zdot)*
           (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
             y*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)),2)/
        ((pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          (pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2)))));
	      
  daei_dxv[5][4] = 
    -((-(((-(xdot*z) + x*zdot)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             (-(ydot*z) + y*zdot)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)))*
           (2*((x*xdot)/mu - (x*xdot + y*ydot + z*zdot)/mu)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             2*((2*xdot*y)/mu - (x*ydot)/mu)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             2*((2*xdot*z)/mu - (x*zdot)/mu)*
              (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu))))/
        (2.*sqrt(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          pow(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                 mu) + x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2),1.5)) + 
       (((2*xdot*y)/mu - (x*ydot)/mu)*(-(ydot*z) + y*zdot) + 
          (-(xdot*z) + x*zdot)*
           ((x*xdot)/mu - (x*xdot + y*ydot + z*zdot)/mu) - 
          z*(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
             x*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)))/
        (sqrt(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2))) - 
       (z*(xdot*z - x*zdot)*
          ((-(xdot*z) + x*zdot)*
             (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu)) + 
            (-(ydot*z) + y*zdot)*
             (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu))))/
        (pow(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2),1.5)*
          sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2))))/
     sqrt(1 - pow((-(xdot*z) + x*zdot)*
           (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
             x*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)) + 
          (-(ydot*z) + y*zdot)*
           (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
             y*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)),2)/
        ((pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          (pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2)))));
	      
  daei_dxv[5][5] = 
    -((-(((-(xdot*z) + x*zdot)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             (-(ydot*z) + y*zdot)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)))*
           (2*(-((xdot*y)/mu) + (2*x*ydot)/mu)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             2*((y*ydot)/mu - (x*xdot + y*ydot + z*zdot)/mu)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             2*((2*ydot*z)/mu - (y*zdot)/mu)*
              (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu))))/
        (2.*sqrt(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          pow(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                 mu) + x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2),1.5)) + 
       ((-((xdot*y)/mu) + (2*x*ydot)/mu)*
           (-(xdot*z) + x*zdot) + 
          (-(ydot*z) + y*zdot)*
           ((y*ydot)/mu - (x*xdot + y*ydot + z*zdot)/mu) - 
          z*(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
             y*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)))/
        (sqrt(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2))) + 
       (z*(-(ydot*z) + y*zdot)*
          ((-(xdot*z) + x*zdot)*
             (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu)) + 
            (-(ydot*z) + y*zdot)*
             (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu))))/
        (pow(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2),1.5)*
          sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2))))/
     sqrt(1 - pow((-(xdot*z) + x*zdot)*
           (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
             x*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)) + 
          (-(ydot*z) + y*zdot)*
           (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
             y*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)),2)/
        ((pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          (pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2)))));
	      
  daei_dxv[5][6] = 
    -((-(((-(xdot*z) + x*zdot)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             (-(ydot*z) + y*zdot)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)))*
           (2*(-((xdot*z)/mu) + (2*x*zdot)/mu)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             2*(-((ydot*z)/mu) + (2*y*zdot)/mu)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) + 
             2*((z*zdot)/mu - (x*xdot + y*ydot + z*zdot)/mu)*
              (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu))))/
        (2.*sqrt(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          pow(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                 mu) + x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2),1.5)) + 
       ((-(xdot*z) + x*zdot)*
           (-((xdot*z)/mu) + (2*x*zdot)/mu) + 
          (-(ydot*z) + y*zdot)*
           (-((ydot*z)/mu) + (2*y*zdot)/mu) + 
          x*(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
             x*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)) + 
          y*(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
             y*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)))/
        (sqrt(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2))) - 
       ((-2*x*(xdot*z - x*zdot) + 2*y*(-(ydot*z) + y*zdot))*
          ((-(xdot*z) + x*zdot)*
             (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu)) + 
            (-(ydot*z) + y*zdot)*
             (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu))))/
        (2.*pow(pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2),1.5)*
          sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2))))/
     sqrt(1 - pow((-(xdot*z) + x*zdot)*
           (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
             x*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)) + 
          (-(ydot*z) + y*zdot)*
           (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
             y*(-(1/
                   sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)),2)/
        ((pow(xdot*z - x*zdot,2) + 
            pow(-(ydot*z) + y*zdot,2))*
          (pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              x*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
              y*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2) + 
            pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
              z*(-(1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
                   + (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu),2)))));
	      
	
/* partials of T (time from periapse) now */

  daei_dxv[6][1] = 
    -((((x*xdot + y*ydot + z*zdot)*
           sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             )*(-2*(-(pow(xdot,2)/mu) + 
                pow(x,2)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - 1/
                 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) - 
             2*(-((xdot*ydot)/mu) + 
                (x*y)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5))*(-((ydot*(x*xdot + y*ydot + z*zdot))/
                   mu) + y*
                 (-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) - 
             2*((x*z)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - (xdot*zdot)/mu)*
              (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu))))/
         (2.*mu*pow(1 - 
             pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2),1.5)) + 
        (2*x*(x*xdot + y*ydot + z*zdot)*
           sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2)))/
         (mu*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*
           sqrt(1 - pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                  mu) + x*
                (-(1/
                     sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                  (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) - 
        ((x*xdot + y*ydot + z*zdot)*
           (2*ydot*(-(xdot*y) + x*ydot) - 
             2*zdot*(xdot*z - x*zdot))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ))/
         (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           sqrt(1 - pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                  mu) + x*
                (-(1/
                     sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                  (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) - 
        (xdot*sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ))/
         (mu*sqrt(1 - pow(-((xdot*
                    (x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) + 
        (-((x*xdot + y*ydot + z*zdot)*
               sqrt(pow(-(xdot*y) + x*ydot,2) + 
                 pow(xdot*z - x*zdot,2) + 
                 pow(-(ydot*z) + y*zdot,2))*
               (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                 (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu)*
               (2*(-(pow(xdot,2)/mu) + 
                    pow(x,2)/
                     pow(pow(x,2) + pow(y,2) + 
                      pow(z,2),1.5) - 
                    1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                      + (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)*
                  (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) + 
                 2*(-((xdot*ydot)/mu) + 
                    (x*y)/
                     pow(pow(x,2) + pow(y,2) + 
                      pow(z,2),1.5))*
                  (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) + 
                 2*((x*z)/
                     pow(pow(x,2) + pow(y,2) + 
                      pow(z,2),1.5) - (xdot*zdot)/mu)*
                  (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu))))/
            (2.*mu*sqrt(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              pow(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2),1.5)) - 
           ((x*xdot + y*ydot + z*zdot)*
              sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)*
              (-2*(-(pow(xdot,2)/mu) + 
                   pow(x,2)/
                    pow(pow(x,2) + pow(y,2) + pow(z,2),
                     1.5) - 
                   1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                    + (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)*
                 (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) - 
                2*(-((xdot*ydot)/mu) + 
                   (x*y)/
                    pow(pow(x,2) + pow(y,2) + pow(z,2),
                     1.5))*
                 (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) - 
                2*((x*z)/
                    pow(pow(x,2) + pow(y,2) + pow(z,2),
                     1.5) - (xdot*zdot)/mu)*
                 (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu))))/
            (2.*mu*pow(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2),1.5)*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) - 
           (2*x*(x*xdot + y*ydot + z*zdot)*
              sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2)))/
            (mu*pow(pow(x,2) + pow(y,2) + pow(z,2),
               1.5)*sqrt(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) + 
           ((x*xdot + y*ydot + z*zdot)*
              (2*ydot*(-(xdot*y) + x*ydot) - 
                2*zdot*(xdot*z - x*zdot))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu))/
            (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              sqrt(1 - pow(-((xdot*
                      (x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) + 
           (xdot*sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu))/
            (mu*sqrt(1 - pow(-((xdot*
                       (x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))))/
         sqrt(1 - (pow(x*xdot + y*ydot + z*zdot,2)*
              (pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              pow(2/
                 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu,2))/
            (pow(mu,2)*(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              (pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2)))))/
      sqrt(mu*pow(2/
           sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
          (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,3)
        )) - (3*mu*x*pow(2/
         sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
        (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,2)*
      (-(((x*xdot + y*ydot + z*zdot)*
             sqrt(pow(-(xdot*y) + x*ydot,2) + 
               pow(xdot*z - x*zdot,2) + 
               pow(-(ydot*z) + y*zdot,2))*
             (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
               (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/
                mu))/
           (mu*sqrt(1 - pow(-((xdot*
                      (x*xdot + y*ydot + z*zdot))/mu) + 
                 x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
               pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                 y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
               pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                 z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2)))) + 
        asin(((x*xdot + y*ydot + z*zdot)*
            sqrt(pow(-(xdot*y) + x*ydot,2) + 
              pow(xdot*z - x*zdot,2) + 
              pow(-(ydot*z) + y*zdot,2))*
            (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
              (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/
               mu))/
          (mu*sqrt(1 - pow(-((xdot*
                     (x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) - 
              pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) - 
              pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2))*
            sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                   mu) + x*
                 (-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) + 
              pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) + 
              pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2))))))/
    (pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*
      pow(mu*pow(2/
           sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
          (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,3)
        ,1.5));
       
  daei_dxv[6][2] = 
    -((((x*xdot + y*ydot + z*zdot)*
           sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             )*(-2*(-((xdot*ydot)/mu) + 
                (x*y)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5))*(-((xdot*(x*xdot + y*ydot + z*zdot))/
                   mu) + x*
                 (-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) - 
             2*(-(pow(ydot,2)/mu) + 
                pow(y,2)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - 1/
                 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) - 
             2*((y*z)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - (ydot*zdot)/mu)*
              (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu))))/
         (2.*mu*pow(1 - 
             pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2),1.5)) + 
        (2*y*(x*xdot + y*ydot + z*zdot)*
           sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2)))/
         (mu*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*
           sqrt(1 - pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                  mu) + x*
                (-(1/
                     sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                  (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) - 
        ((x*xdot + y*ydot + z*zdot)*
           (-2*xdot*(-(xdot*y) + x*ydot) + 
             2*zdot*(-(ydot*z) + y*zdot))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ))/
         (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           sqrt(1 - pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                  mu) + x*
                (-(1/
                     sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                  (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) - 
        (ydot*sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ))/
         (mu*sqrt(1 - pow(-((xdot*
                    (x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) + 
        (-((x*xdot + y*ydot + z*zdot)*
               sqrt(pow(-(xdot*y) + x*ydot,2) + 
                 pow(xdot*z - x*zdot,2) + 
                 pow(-(ydot*z) + y*zdot,2))*
               (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                 (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu)*
               (2*(-((xdot*ydot)/mu) + 
                    (x*y)/
                     pow(pow(x,2) + pow(y,2) + 
                      pow(z,2),1.5))*
                  (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) + 
                 2*(-(pow(ydot,2)/mu) + 
                    pow(y,2)/
                     pow(pow(x,2) + pow(y,2) + 
                      pow(z,2),1.5) - 
                    1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                      + (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)*
                  (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) + 
                 2*((y*z)/
                     pow(pow(x,2) + pow(y,2) + 
                      pow(z,2),1.5) - (ydot*zdot)/mu)*
                  (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu))))/
            (2.*mu*sqrt(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              pow(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2),1.5)) - 
           ((x*xdot + y*ydot + z*zdot)*
              sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)*
              (-2*(-((xdot*ydot)/mu) + 
                   (x*y)/
                    pow(pow(x,2) + pow(y,2) + pow(z,2),
                     1.5))*
                 (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) - 
                2*(-(pow(ydot,2)/mu) + 
                   pow(y,2)/
                    pow(pow(x,2) + pow(y,2) + pow(z,2),
                     1.5) - 
                   1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                    + (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)*
                 (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) - 
                2*((y*z)/
                    pow(pow(x,2) + pow(y,2) + pow(z,2),
                     1.5) - (ydot*zdot)/mu)*
                 (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu))))/
            (2.*mu*pow(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2),1.5)*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) - 
           (2*y*(x*xdot + y*ydot + z*zdot)*
              sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2)))/
            (mu*pow(pow(x,2) + pow(y,2) + pow(z,2),
               1.5)*sqrt(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) + 
           ((x*xdot + y*ydot + z*zdot)*
              (-2*xdot*(-(xdot*y) + x*ydot) + 
                2*zdot*(-(ydot*z) + y*zdot))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu))/
            (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              sqrt(1 - pow(-((xdot*
                      (x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) + 
           (ydot*sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu))/
            (mu*sqrt(1 - pow(-((xdot*
                       (x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))))/
         sqrt(1 - (pow(x*xdot + y*ydot + z*zdot,2)*
              (pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              pow(2/
                 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu,2))/
            (pow(mu,2)*(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              (pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2)))))/
      sqrt(mu*pow(2/
           sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
          (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,3)
        )) - (3*mu*y*pow(2/
         sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
        (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,2)*
      (-(((x*xdot + y*ydot + z*zdot)*
             sqrt(pow(-(xdot*y) + x*ydot,2) + 
               pow(xdot*z - x*zdot,2) + 
               pow(-(ydot*z) + y*zdot,2))*
             (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
               (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/
                mu))/
           (mu*sqrt(1 - pow(-((xdot*
                      (x*xdot + y*ydot + z*zdot))/mu) + 
                 x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
               pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                 y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
               pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                 z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2)))) + 
        asin(((x*xdot + y*ydot + z*zdot)*
            sqrt(pow(-(xdot*y) + x*ydot,2) + 
              pow(xdot*z - x*zdot,2) + 
              pow(-(ydot*z) + y*zdot,2))*
            (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
              (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/
               mu))/
          (mu*sqrt(1 - pow(-((xdot*
                     (x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) - 
              pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) - 
              pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2))*
            sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                   mu) + x*
                 (-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) + 
              pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) + 
              pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2))))))/
    (pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*
      pow(mu*pow(2/
           sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
          (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,3)
        ,1.5));

  daei_dxv[6][3] = 
    -((((x*xdot + y*ydot + z*zdot)*
           sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             )*(-2*((x*z)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - (xdot*zdot)/mu)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) - 
             2*((y*z)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - (ydot*zdot)/mu)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) - 
             2*(pow(z,2)/
                 pow(pow(x,2) + pow(y,2) + pow(z,2),
                  1.5) - 1/
                 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                pow(zdot,2)/mu + 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)*
              (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu))))/
         (2.*mu*pow(1 - 
             pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2),1.5)) + 
        (2*z*(x*xdot + y*ydot + z*zdot)*
           sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2)))/
         (mu*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*
           sqrt(1 - pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                  mu) + x*
                (-(1/
                     sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                  (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) - 
        ((x*xdot + y*ydot + z*zdot)*
           (2*xdot*(xdot*z - x*zdot) - 
             2*ydot*(-(ydot*z) + y*zdot))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ))/
         (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           sqrt(1 - pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                  mu) + x*
                (-(1/
                     sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                  (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) - 
        (zdot*sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ))/
         (mu*sqrt(1 - pow(-((xdot*
                    (x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) + 
        (-((x*xdot + y*ydot + z*zdot)*
               sqrt(pow(-(xdot*y) + x*ydot,2) + 
                 pow(xdot*z - x*zdot,2) + 
                 pow(-(ydot*z) + y*zdot,2))*
               (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                 (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu)*
               (2*((x*z)/
                     pow(pow(x,2) + pow(y,2) + 
                      pow(z,2),1.5) - (xdot*zdot)/mu)*
                  (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) + 
                 2*((y*z)/
                     pow(pow(x,2) + pow(y,2) + 
                      pow(z,2),1.5) - (ydot*zdot)/mu)*
                  (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) + 
                 2*(pow(z,2)/
                     pow(pow(x,2) + pow(y,2) + 
                      pow(z,2),1.5) - 
                    1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                      - pow(zdot,2)/mu + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)*
                  (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu))))/
            (2.*mu*sqrt(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              pow(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2),1.5)) - 
           ((x*xdot + y*ydot + z*zdot)*
              sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)*
              (-2*((x*z)/
                    pow(pow(x,2) + pow(y,2) + pow(z,2),
                     1.5) - (xdot*zdot)/mu)*
                 (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) - 
                2*((y*z)/
                    pow(pow(x,2) + pow(y,2) + pow(z,2),
                     1.5) - (ydot*zdot)/mu)*
                 (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) - 
                2*(pow(z,2)/
                    pow(pow(x,2) + pow(y,2) + pow(z,2),
                     1.5) - 
                   1/
                    sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                    - pow(zdot,2)/mu + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)*
                 (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu))))/
            (2.*mu*pow(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2),1.5)*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) - 
           (2*z*(x*xdot + y*ydot + z*zdot)*
              sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2)))/
            (mu*pow(pow(x,2) + pow(y,2) + pow(z,2),
               1.5)*sqrt(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) + 
           ((x*xdot + y*ydot + z*zdot)*
              (2*xdot*(xdot*z - x*zdot) - 
                2*ydot*(-(ydot*z) + y*zdot))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu))/
            (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              sqrt(1 - pow(-((xdot*
                      (x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) + 
           (zdot*sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu))/
            (mu*sqrt(1 - pow(-((xdot*
                       (x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))))/
         sqrt(1 - (pow(x*xdot + y*ydot + z*zdot,2)*
              (pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              pow(2/
                 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu,2))/
            (pow(mu,2)*(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              (pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2)))))/
      sqrt(mu*pow(2/
           sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
          (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,3)
        )) - (3*mu*z*pow(2/
         sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
        (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,2)*
      (-(((x*xdot + y*ydot + z*zdot)*
             sqrt(pow(-(xdot*y) + x*ydot,2) + 
               pow(xdot*z - x*zdot,2) + 
               pow(-(ydot*z) + y*zdot,2))*
             (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
               (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/
                mu))/
           (mu*sqrt(1 - pow(-((xdot*
                      (x*xdot + y*ydot + z*zdot))/mu) + 
                 x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
               pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                 y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
               pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                 z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2)))) + 
        asin(((x*xdot + y*ydot + z*zdot)*
            sqrt(pow(-(xdot*y) + x*ydot,2) + 
              pow(xdot*z - x*zdot,2) + 
              pow(-(ydot*z) + y*zdot,2))*
            (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
              (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/
               mu))/
          (mu*sqrt(1 - pow(-((xdot*
                     (x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) - 
              pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) - 
              pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2))*
            sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                   mu) + x*
                 (-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) + 
              pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) + 
              pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2))))))/
    (pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*
      pow(mu*pow(2/
           sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
          (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,3)
        ,1.5));
       
  daei_dxv[6][4] = 
    -((((x*xdot + y*ydot + z*zdot)*
           sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             )*(-2*((x*xdot)/mu - 
                (x*xdot + y*ydot + z*zdot)/mu)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) - 
             2*((2*xdot*y)/mu - (x*ydot)/mu)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) - 
             2*((2*xdot*z)/mu - (x*zdot)/mu)*
              (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu))))/
         (2.*mu*pow(1 - 
             pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2),1.5)) + 
        (2*xdot*(x*xdot + y*ydot + z*zdot)*
           sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2)))/
         (pow(mu,2)*sqrt(1 - 
             pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) - 
        ((x*xdot + y*ydot + z*zdot)*
           (-2*y*(-(xdot*y) + x*ydot) + 2*z*(xdot*z - x*zdot))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ))/
         (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           sqrt(1 - pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                  mu) + x*
                (-(1/
                     sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                  (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) - 
        (x*sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ))/
         (mu*sqrt(1 - pow(-((xdot*
                    (x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) + 
        (-((x*xdot + y*ydot + z*zdot)*
               sqrt(pow(-(xdot*y) + x*ydot,2) + 
                 pow(xdot*z - x*zdot,2) + 
                 pow(-(ydot*z) + y*zdot,2))*
               (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                 (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu)*
               (2*((x*xdot)/mu - 
                    (x*xdot + y*ydot + z*zdot)/mu)*
                  (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) + 
                 2*((2*xdot*y)/mu - (x*ydot)/mu)*
                  (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) + 
                 2*((2*xdot*z)/mu - (x*zdot)/mu)*
                  (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu))))/
            (2.*mu*sqrt(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              pow(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2),1.5)) - 
           ((x*xdot + y*ydot + z*zdot)*
              sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)*
              (-2*((x*xdot)/mu - 
                   (x*xdot + y*ydot + z*zdot)/mu)*
                 (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) - 
                2*((2*xdot*y)/mu - (x*ydot)/mu)*
                 (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) - 
                2*((2*xdot*z)/mu - (x*zdot)/mu)*
                 (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu))))/
            (2.*mu*pow(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2),1.5)*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) - 
           (2*xdot*(x*xdot + y*ydot + z*zdot)*
              sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2)))/
            (pow(mu,2)*sqrt(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) + 
           ((x*xdot + y*ydot + z*zdot)*
              (-2*y*(-(xdot*y) + x*ydot) + 
                2*z*(xdot*z - x*zdot))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu))/
            (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              sqrt(1 - pow(-((xdot*
                      (x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) + 
           (x*sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu))/
            (mu*sqrt(1 - pow(-((xdot*
                       (x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))))/
         sqrt(1 - (pow(x*xdot + y*ydot + z*zdot,2)*
              (pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              pow(2/
                 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu,2))/
            (pow(mu,2)*(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              (pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2)))))/
      sqrt(mu*pow(2/
           sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
          (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,3)
        )) - (3*xdot*pow(2/
         sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
        (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,2)*
      (-(((x*xdot + y*ydot + z*zdot)*
             sqrt(pow(-(xdot*y) + x*ydot,2) + 
               pow(xdot*z - x*zdot,2) + 
               pow(-(ydot*z) + y*zdot,2))*
             (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
               (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/
                mu))/
           (mu*sqrt(1 - pow(-((xdot*
                      (x*xdot + y*ydot + z*zdot))/mu) + 
                 x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
               pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                 y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
               pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                 z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2)))) + 
        asin(((x*xdot + y*ydot + z*zdot)*
            sqrt(pow(-(xdot*y) + x*ydot,2) + 
              pow(xdot*z - x*zdot,2) + 
              pow(-(ydot*z) + y*zdot,2))*
            (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
              (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/
               mu))/
          (mu*sqrt(1 - pow(-((xdot*
                     (x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) - 
              pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) - 
              pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2))*
            sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                   mu) + x*
                 (-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) + 
              pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) + 
              pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2))))))/
    pow(mu*pow(2/
         sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
        (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,3),
     1.5);
     
       
  daei_dxv[6][5] = 
    -((((x*xdot + y*ydot + z*zdot)*
           sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             )*(-2*(-((xdot*y)/mu) + (2*x*ydot)/mu)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) - 
             2*((y*ydot)/mu - (x*xdot + y*ydot + z*zdot)/mu)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) - 
             2*((2*ydot*z)/mu - (y*zdot)/mu)*
              (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu))))/
         (2.*mu*pow(1 - 
             pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2),1.5)) + 
        (2*ydot*(x*xdot + y*ydot + z*zdot)*
           sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2)))/
         (pow(mu,2)*sqrt(1 - 
             pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) - 
        ((x*xdot + y*ydot + z*zdot)*
           (2*x*(-(xdot*y) + x*ydot) - 
             2*z*(-(ydot*z) + y*zdot))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ))/
         (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           sqrt(1 - pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                  mu) + x*
                (-(1/
                     sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                  (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) - 
        (y*sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ))/
         (mu*sqrt(1 - pow(-((xdot*
                    (x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) + 
        (-((x*xdot + y*ydot + z*zdot)*
               sqrt(pow(-(xdot*y) + x*ydot,2) + 
                 pow(xdot*z - x*zdot,2) + 
                 pow(-(ydot*z) + y*zdot,2))*
               (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                 (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu)*
               (2*(-((xdot*y)/mu) + (2*x*ydot)/mu)*
                  (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) + 
                 2*((y*ydot)/mu - 
                    (x*xdot + y*ydot + z*zdot)/mu)*
                  (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) + 
                 2*((2*ydot*z)/mu - (y*zdot)/mu)*
                  (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu))))/
            (2.*mu*sqrt(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              pow(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2),1.5)) - 
           ((x*xdot + y*ydot + z*zdot)*
              sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)*
              (-2*(-((xdot*y)/mu) + (2*x*ydot)/mu)*
                 (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) - 
                2*((y*ydot)/mu - 
                   (x*xdot + y*ydot + z*zdot)/mu)*
                 (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) - 
                2*((2*ydot*z)/mu - (y*zdot)/mu)*
                 (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu))))/
            (2.*mu*pow(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2),1.5)*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) - 
           (2*ydot*(x*xdot + y*ydot + z*zdot)*
              sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2)))/
            (pow(mu,2)*sqrt(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) + 
           ((x*xdot + y*ydot + z*zdot)*
              (2*x*(-(xdot*y) + x*ydot) - 
                2*z*(-(ydot*z) + y*zdot))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu))/
            (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              sqrt(1 - pow(-((xdot*
                      (x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) + 
           (y*sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu))/
            (mu*sqrt(1 - pow(-((xdot*
                       (x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))))/
         sqrt(1 - (pow(x*xdot + y*ydot + z*zdot,2)*
              (pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              pow(2/
                 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu,2))/
            (pow(mu,2)*(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              (pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2)))))/
      sqrt(mu*pow(2/
           sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
          (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,3)
        )) - (3*ydot*pow(2/
         sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
        (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,2)*
      (-(((x*xdot + y*ydot + z*zdot)*
             sqrt(pow(-(xdot*y) + x*ydot,2) + 
               pow(xdot*z - x*zdot,2) + 
               pow(-(ydot*z) + y*zdot,2))*
             (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
               (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/
                mu))/
           (mu*sqrt(1 - pow(-((xdot*
                      (x*xdot + y*ydot + z*zdot))/mu) + 
                 x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
               pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                 y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
               pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                 z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2)))) + 
        asin(((x*xdot + y*ydot + z*zdot)*
            sqrt(pow(-(xdot*y) + x*ydot,2) + 
              pow(xdot*z - x*zdot,2) + 
              pow(-(ydot*z) + y*zdot,2))*
            (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
              (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/
               mu))/
          (mu*sqrt(1 - pow(-((xdot*
                     (x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) - 
              pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) - 
              pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2))*
            sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                   mu) + x*
                 (-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) + 
              pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) + 
              pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2))))))/
    pow(mu*pow(2/
         sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
        (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,3),
     1.5);
     
  daei_dxv[6][6] = 
    -((((x*xdot + y*ydot + z*zdot)*
           sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             )*(-2*(-((xdot*z)/mu) + (2*x*zdot)/mu)*
              (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) - 
             2*(-((ydot*z)/mu) + (2*y*zdot)/mu)*
              (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu)) - 
             2*((z*zdot)/mu - (x*xdot + y*ydot + z*zdot)/mu)*
              (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu))))/
         (2.*mu*pow(1 - 
             pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2),1.5)) + 
        (2*zdot*(x*xdot + y*ydot + z*zdot)*
           sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2)))/
         (pow(mu,2)*sqrt(1 - 
             pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) - 
        ((x*xdot + y*ydot + z*zdot)*
           (-2*x*(xdot*z - x*zdot) + 2*y*(-(ydot*z) + y*zdot))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ))/
         (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           sqrt(1 - pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                  mu) + x*
                (-(1/
                     sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                  (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) - 
        (z*sqrt(pow(-(xdot*y) + x*ydot,2) + 
             pow(xdot*z - x*zdot,2) + 
             pow(-(ydot*z) + y*zdot,2))*
           (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
             (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu
             ))/
         (mu*sqrt(1 - pow(-((xdot*
                    (x*xdot + y*ydot + z*zdot))/mu) + 
               x*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
               y*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2) - 
             pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
               z*(-(1/
                     sqrt(pow(x,2) + pow(y,2) + pow(z,2))
                     ) + (pow(xdot,2) + pow(ydot,2) + 
                     pow(zdot,2))/mu),2))) + 
        (-((x*xdot + y*ydot + z*zdot)*
               sqrt(pow(-(xdot*y) + x*ydot,2) + 
                 pow(xdot*z - x*zdot,2) + 
                 pow(-(ydot*z) + y*zdot,2))*
               (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                 (pow(xdot,2) + pow(ydot,2) + 
                    pow(zdot,2))/mu)*
               (2*(-((xdot*z)/mu) + (2*x*zdot)/mu)*
                  (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) + 
                 2*(-((ydot*z)/mu) + (2*y*zdot)/mu)*
                  (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) + 
                 2*((z*zdot)/mu - 
                    (x*xdot + y*ydot + z*zdot)/mu)*
                  (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                    z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                       (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu))))/
            (2.*mu*sqrt(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              pow(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2),1.5)) - 
           ((x*xdot + y*ydot + z*zdot)*
              sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu)*
              (-2*(-((xdot*z)/mu) + (2*x*zdot)/mu)*
                 (-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) - 
                2*(-((ydot*z)/mu) + (2*y*zdot)/mu)*
                 (-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu)) - 
                2*((z*zdot)/mu - 
                   (x*xdot + y*ydot + z*zdot)/mu)*
                 (-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                   z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                      (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu))))/
            (2.*mu*pow(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2),1.5)*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) - 
           (2*zdot*(x*xdot + y*ydot + z*zdot)*
              sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2)))/
            (pow(mu,2)*sqrt(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) + 
           ((x*xdot + y*ydot + z*zdot)*
              (-2*x*(xdot*z - x*zdot) + 
                2*y*(-(ydot*z) + y*zdot))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu))/
            (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              sqrt(1 - pow(-((xdot*
                      (x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))) + 
           (z*sqrt(pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu))/
            (mu*sqrt(1 - pow(-((xdot*
                       (x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))))/
         sqrt(1 - (pow(x*xdot + y*ydot + z*zdot,2)*
              (pow(-(xdot*y) + x*ydot,2) + 
                pow(xdot*z - x*zdot,2) + 
                pow(-(ydot*z) + y*zdot,2))*
              pow(2/
                 sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
                (pow(xdot,2) + pow(ydot,2) + 
                   pow(zdot,2))/mu,2))/
            (pow(mu,2)*(1 - 
                pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2))*
              (pow(-((xdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                  x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((ydot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) + 
                pow(-((zdot*(x*xdot + y*ydot + z*zdot))/
                     mu) + 
                  z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                     (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2)))))/
      sqrt(mu*pow(2/
           sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
          (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,3)
        )) - (3*zdot*pow(2/
         sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
        (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,2)*
      (-(((x*xdot + y*ydot + z*zdot)*
             sqrt(pow(-(xdot*y) + x*ydot,2) + 
               pow(xdot*z - x*zdot,2) + 
               pow(-(ydot*z) + y*zdot,2))*
             (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
               (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/
                mu))/
           (mu*sqrt(1 - pow(-((xdot*
                      (x*xdot + y*ydot + z*zdot))/mu) + 
                 x*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
               pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                 y*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2) - 
               pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                 z*(-(1/
                       sqrt(pow(x,2) + pow(y,2) + 
                       pow(z,2))) + 
                    (pow(xdot,2) + pow(ydot,2) + 
                       pow(zdot,2))/mu),2)))) + 
        asin(((x*xdot + y*ydot + z*zdot)*
            sqrt(pow(-(xdot*y) + x*ydot,2) + 
              pow(xdot*z - x*zdot,2) + 
              pow(-(ydot*z) + y*zdot,2))*
            (2/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
              (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/
               mu))/
          (mu*sqrt(1 - pow(-((xdot*
                     (x*xdot + y*ydot + z*zdot))/mu) + 
                x*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) - 
              pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) - 
              pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2))*
            sqrt(pow(-((xdot*(x*xdot + y*ydot + z*zdot))/
                   mu) + x*
                 (-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) + 
              pow(-((ydot*(x*xdot + y*ydot + z*zdot))/mu) + 
                y*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2) + 
              pow(-((zdot*(x*xdot + y*ydot + z*zdot))/mu) + 
                z*(-(1/
                      sqrt(pow(x,2) + pow(y,2) + 
                      pow(z,2))) + 
                   (pow(xdot,2) + pow(ydot,2) + 
                      pow(zdot,2))/mu),2))))))/
    pow(mu*pow(2/
         sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - 
        (pow(xdot,2) + pow(ydot,2) + pow(zdot,2))/mu,3),
     1.5);
     
  return;
}
  
                                    	      	      	      	      	      	      



















































                                         	  	  	  	  	  	  

	 	 	 	 	 	                             
       
                     
