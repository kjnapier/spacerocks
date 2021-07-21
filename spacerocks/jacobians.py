import numpy as np


def xyz_to_ecl_jacobian(x, y, z):

    dλ_dx = -y/(x**2 + y**2)
    dλ_dy = x/(x**2 + y**2)

    dβ_dx = - np.sign(x) * np.sign(z) * np.sqrt(((x**2 * z**2) / ((x**2 + y**2)*(x**2 + y**2 + z**2)**2)))
    dβ_dy = - np.sign(y) * np.sign(z) * np.sqrt(((y**2 * z**2) / ((x**2 + y**2)*(x**2 + y**2 + z**2)**2)))
    dβ_dz = np.sqrt(((x**2 + y**2) / (x**2 + y**2 + z**2)**2))

    jacobian = np.array([[dλ_dx, dλ_dy, 0, 0, 0, 0],
                         [dβ_dx, dβ_dy, dβ_dz, 0, 0, 0]])

    return jacobian


def ecl_to_equa_jacobian(ra, dec, epsilon):

    dβ_ddec = (np.cos(dec)*np.cos(epsilon) + np.sin(ra)*np.sin(dec)*np.sin(epsilon))**2 / \
              (1 - (np.cos(epsilon)*np.sin(dec) - np.cos(dec)*np.sin(ra)*np.sin(epsilon))**2)

    dβ_dra = (np.cos(ra) * np.cos(dec) * np.sin(epsilon))**2 / \
             (1 - (np.cos(epsilon)*np.sin(dec) - np.cos(dec)*np.sin(ra)*np.sin(epsilon))**2)

    dλ_ddec = (np.cos(ra)**2 * np.sec(dec)**4 * np.sin(epsilon)**2) \
              / (np.cos(ra)**2 + (np.cos(epsilon)*np.sin(ra) + np.sin(epsilon)*np.tan(dec))**2)**2

    dλ_dra = (np.cos(epsilon) + np.sin(ra) * np.sin(epsilon) * np.tan(dec))**2 \
             / ((np.cos(ra)**2 + np.cos(epsilon)**2 * np.sin(ra)**2) \
             + np.tan(dec) * (np.sin(ra) * np.sin(2*epsilon) \
             + np.sin(epsilon)**2 * np.tan(dec))**2)

    return np.array([[dλ_dra, dλ_ddec], [dβ_dra, dβ_ddec]])


def kep_to_xyz_jacobian(x, y, z, vx, vy, vz, mu):

    jacobian = np.zeros([6, 6])

    jacobian[0, 0] = \
       (2*x)/(pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*\
	   pow(2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
	       (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,2));

    jacobian[0, 1] = \
        (2*y)/(pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*\
	    pow(2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
		(pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,2));

    jacobian[0, 2] = \
        (2*z)/(pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*\
	    pow(2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
		(pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,2));

    jacobian[0, 3] = \
            (2*vx)/(mu*pow(2/\
			 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
			 (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,2));

    jacobian[0, 4] = \
            (2*vy)/(mu*pow(2/\
			 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
			 (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,2));

    jacobian[0, 5] = \
            (2*vz)/(mu*pow(2/\
		     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
		     (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,2));

    # Partials for e now */

    jacobian[1, 0] = \
    (2*(-(pow(vx,2)/mu) + \
        pow(x,2)/\
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - \
        1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + \
        (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)*\
      (-((vx*(x*vx + y*vy + z*vz))/mu) + \
        x*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
       + 2*(-((vx*vy)/mu) + \
        (x*y)/pow(pow(x,2) + pow(y,2) + pow(z,2),1.5))*\
      (-((vy*(x*vx + y*vy + z*vz))/mu) + \
        y*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
       + 2*((x*z)/\
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - \
        (vx*vz)/mu)*(-((vz*(x*vx + y*vy + z*vz))/\
           mu) + z*(-(1/\
              np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
     )/(2.*np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
            mu) + x*(-(1/\
               np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2) + pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
         y*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2) + pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
         z*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2)));

    jacobian[1, 1] = \
    (2*(-((vx*vy)/mu) + (x*y)/\
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5))*\
      (-((vx*(x*vx + y*vy + z*vz))/mu) + \
        x*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
       + 2*(-(pow(vy,2)/mu) + \
        pow(y,2)/\
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - \
        1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + \
        (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)*\
      (-((vy*(x*vx + y*vy + z*vz))/mu) + \
        y*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
       + 2*((y*z)/\
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - \
        (vy*vz)/mu)*(-((vz*(x*vx + y*vy + z*vz))/\
           mu) + z*(-(1/\
              np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
     )/(2.*np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
            mu) + x*(-(1/\
               np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2) + pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
         y*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2) + pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
         z*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2)));

    jacobian[1, 2] = \
    (2*((x*z)/pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - \
        (vx*vz)/mu)*(-((vx*(x*vx + y*vy + z*vz))/\
           mu) + x*(-(1/\
              np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
       + 2*((y*z)/\
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - \
        (vy*vz)/mu)*(-((vy*(x*vx + y*vy + z*vz))/\
           mu) + y*(-(1/\
              np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
       + 2*(pow(z,2)/\
         pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - \
        1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
        pow(vz,2)/mu + \
        (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)*\
      (-((vz*(x*vx + y*vy + z*vz))/mu) + \
        z*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
     )/(2.*np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
            mu) + x*(-(1/\
               np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2) + pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
         y*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2) + pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
         z*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2)));

    jacobian[1, 3] = \
    (2*((x*vx)/mu - (x*vx + y*vy + z*vz)/mu)*\
      (-((vx*(x*vx + y*vy + z*vz))/mu) + \
        x*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
       + 2*((2*vx*y)/mu - (x*vy)/mu)*\
      (-((vy*(x*vx + y*vy + z*vz))/mu) + \
        y*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
       + 2*((2*vx*z)/mu - (x*vz)/mu)*\
      (-((vz*(x*vx + y*vy + z*vz))/mu) + \
        z*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
     )/(2.*np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
            mu) + x*(-(1/\
               np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2) + pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
         y*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2) + pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
         z*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2)));

    jacobian[1, 4] = \
    (2*(-((vx*y)/mu) + (2*x*vy)/mu)*\
      (-((vx*(x*vx + y*vy + z*vz))/mu) + \
        x*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
       + 2*((y*vy)/mu - (x*vx + y*vy + z*vz)/mu)*\
      (-((vy*(x*vx + y*vy + z*vz))/mu) + \
        y*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
       + 2*((2*vy*z)/mu - (y*vz)/mu)*\
      (-((vz*(x*vx + y*vy + z*vz))/mu) + \
        z*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
     )/(2.*np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
            mu) + x*(-(1/\
               np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2) + pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
         y*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2) + pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
         z*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2)));

    jacobian[1, 5] = \
    (2*(-((vx*z)/mu) + (2*x*vz)/mu)*\
      (-((vx*(x*vx + y*vy + z*vz))/mu) + \
        x*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
       + 2*(-((vy*z)/mu) + (2*y*vz)/mu)*\
      (-((vy*(x*vx + y*vy + z*vz))/mu) + \
        y*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
       + 2*((z*vz)/mu - (x*vx + y*vy + z*vz)/mu)*\
      (-((vz*(x*vx + y*vy + z*vz))/mu) + \
        z*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
           (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu))\
     )/(2.*np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
            mu) + x*(-(1/\
               np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) +\
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2) + pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
         y*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2) + pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
         z*(-(1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))) + \
            (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu)\
         ,2)));

    # Partials of i now */

    jacobian[2, 0] = \
    -((-((-(vx*y) + x*vy)*\
           (2*vy*(-(vx*y) + x*vy) - \
             2*vz*(vx*z - x*vz)))/\
        (2.*pow(pow(-(vx*y) + x*vy,2) + \
            pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2),1.5)) + \
       vy/np.sqrt(pow(-(vx*y) + x*vy,2) + \
          pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2)\
          ))/\
     np.sqrt(1 - pow(-(vx*y) + x*vy,2)/\
        (pow(-(vx*y) + x*vy,2) + \
          pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2)\
          )));

    jacobian[2, 1] = \
    -((-((-(vx*y) + x*vy)*\
           (-2*vx*(-(vx*y) + x*vy) + \
             2*vz*(-(vy*z) + y*vz)))/\
        (2.*pow(pow(-(vx*y) + x*vy,2) + \
            pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2),1.5)) - \
       vx/np.sqrt(pow(-(vx*y) + x*vy,2) + \
          pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2)\
          ))/\
     np.sqrt(1 - pow(-(vx*y) + x*vy,2)/\
        (pow(-(vx*y) + x*vy,2) + \
          pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2)\
          )));

    jacobian[2, 2] = \
    ((-(vx*y) + x*vy)*(2*vx*(vx*z - x*vz) - \
       2*vy*(-(vy*z) + y*vz)))/\
    (2.*pow(pow(-(vx*y) + x*vy,2) + \
       pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2),\
      1.5)*np.sqrt(1 - pow(-(vx*y) + x*vy,2)/\
        (pow(-(vx*y) + x*vy,2) + \
          pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2)\
          )));

    jacobian[2, 3] = \
    -((-((-(vx*y) + x*vy)*\
           (-2*y*(-(vx*y) + x*vy) + 2*z*(vx*z - x*vz)))/
        (2.*pow(pow(-(vx*y) + x*vy,2) + \
            pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2),1.5)) - \
       y/np.sqrt(pow(-(vx*y) + x*vy,2) + \
          pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2)\
          ))/\
     np.sqrt(1 - pow(-(vx*y) + x*vy,2)/\
        (pow(-(vx*y) + x*vy,2) + \
          pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2)\
          )));

    jacobian[2, 4] = \
    -((-((-(vx*y) + x*vy)*\
           (2*x*(-(vx*y) + x*vy) - 2*z*(-(vy*z) + y*vz))\
           )/\
        (2.*pow(pow(-(vx*y) + x*vy,2) + \
            pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2),1.5)) + \
       x/np.sqrt(pow(-(vx*y) + x*vy,2) + \
          pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2)\
          ))/\
     np.sqrt(1 - pow(-(vx*y) + x*vy,2)/\
        (pow(-(vx*y) + x*vy,2) + \
          pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2)\
          )));

    jacobian[2, 5] = \
    ((-(vx*y) + x*vy)*(-2*x*(vx*z - x*vz) + \
       2*y*(-(vy*z) + y*vz)))/\
     (2.*pow(pow(-(vx*y) + x*vy,2) + \
       pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2),\
      1.5)*np.sqrt(1 - pow(-(vx*y) + x*vy,2)/\
        (pow(-(vx*y) + x*vy,2) + \
          pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2)\
          )));


    # Partials of capital Omega (long of asc node) now */

    jacobian[3, 0] = \
    -(((vz*(vx*z - x*vz)*(-(vx*z) + x*vz))/\
        pow(pow(vx*z - x*vz,2) + \
          pow(-(vy*z) + y*vz,2),1.5) + \
       vz/np.sqrt(pow(vx*z - x*vz,2) + \
          pow(-(vy*z) + y*vz,2)))/\
     np.sqrt(1 - pow(-(vx*z) + x*vz,2)/\
        (pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2))\
       ));

    jacobian[3, 1] = \
    (vz*(-(vx*z) + x*vz)*(-(vy*z) + y*vz))/\
     (pow(pow(vx*z - x*vz,2) + \
       pow(-(vy*z) + y*vz,2),1.5)*\
     np.sqrt(1 - pow(-(vx*z) + x*vz,2)/\
        (pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2))\
       ));

    jacobian[3, 2] = \
    -((-((-(vx*z) + x*vz)*\
           (2*vx*(vx*z - x*vz) - \
             2*vy*(-(vy*z) + y*vz)))/\
        (2.*pow(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2),1.5)) - \
       vx/np.sqrt(pow(vx*z - x*vz,2) + \
          pow(-(vy*z) + y*vz,2)))/\
     np.sqrt(1 - pow(-(vx*z) + x*vz,2)/\
        (pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2))\
       ));

    jacobian[3, 3] = \
    -((-((z*(vx*z - x*vz)*(-(vx*z) + x*vz))/\
          pow(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2),1.5)) - \
       z/np.sqrt(pow(vx*z - x*vz,2) + \
          pow(-(vy*z) + y*vz,2)))/\
     np.sqrt(1 - pow(-(vx*z) + x*vz,2)/\
        (pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2))\
       ));

    jacobian[3, 4] = \
    -((z*(-(vx*z) + x*vz)*(-(vy*z) + y*vz))/\
     (pow(pow(vx*z - x*vz,2) + \
         pow(-(vy*z) + y*vz,2),1.5)*\
       np.sqrt(1 - pow(-(vx*z) + x*vz,2)/\
          (pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2)))));

    jacobian[3, 5] = \
    -((-((-(vx*z) + x*vz)*\
           (-2*x*(vx*z - x*vz) + 2*y*(-(vy*z) + y*vz)))/\
        (2.*pow(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2),1.5)) + \
       x/np.sqrt(pow(vx*z - x*vz,2) + \
          pow(-(vy*z) + y*vz,2)))/\
     np.sqrt(1 - pow(-(vx*z) + x*vz,2)/\
        (pow(vx*z - x*vz,2) + pow(-(vy*z) + y*vz,2))\
       ));


    # partials of small omega (arg of per) now */

    jacobian[4, 0] = \
    -((-(((-(vx*z) + x*vz)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             (-(vy*z) + y*vz)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)))*\
           (2*(-(pow(vx,2)/mu) + \
                pow(x,2)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - 1/\
                 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) +\
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             2*(-((vx*vy)/mu) + \
                (x*y)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5))*(-((vy*(x*vx + y*vy + z*vz))/\
                   mu) + y*\
                 (-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             2*((x*z)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - (vx*vz)/mu)*\
              (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu))))/\
        (2.*np.sqrt(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          pow(pow(-((vx*(x*vx + y*vy + z*vz))/\
                 mu) + x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2),1.5)) + \
       ((-((vx*vy)/mu) + \
             (x*y)/\
              pow(pow(x,2) + pow(y,2) + pow(z,2),1.5))*\
           (-(vy*z) + y*vz) + \
          (-(vx*z) + x*vz)*\
           (-(pow(vx,2)/mu) + \
             pow(x,2)/\
              pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)\
              - 1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ) + vz*(-((vx*(x*vx + y*vy + z*vz))/\
                mu) + x*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)))/\
        (np.sqrt(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2))) + \
       (vz*(vx*z - x*vz)*\
          ((-(vx*z) + x*vz)*\
             (-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu)) + \
            (-(vy*z) + y*vz)*\
             (-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu))))/\
        (pow(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2),1.5)*\
          np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2))))/\
     np.sqrt(1 - pow((-(vx*z) + x*vz)*\
           (-((vx*(x*vx + y*vy + z*vz))/mu) + \
             x*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)) + \
          (-(vy*z) + y*vz)*\
           (-((vy*(x*vx + y*vy + z*vz))/mu) + \
             y*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)),2)/\
        ((pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          (pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2)))));

    jacobian[4, 1] = \
    -((-(((-(vx*z) + x*vz)*
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/
                      np.sqrt(pow(x,2) + pow(y,2) +\
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             (-(vy*z) + y*vz)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)))*\
           (2*(-((vx*vy)/mu) + \
                (x*y)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5))*(-((vx*(x*vx + y*vy + z*vz))/\
                   mu) + x*\
                 (-(1/
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             2*(-(pow(vy,2)/mu) + \
                pow(y,2)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - 1/\
                 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             2*((y*z)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - (vy*vz)/mu)*\
              (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu))))/\
        (2.*np.sqrt(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          pow(pow(-((vx*(x*vx + y*vy + z*vz))/\
                 mu) + x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2),1.5)) + \
       ((-((vx*vy)/mu) + \
             (x*y)/\
              pow(pow(x,2) + pow(y,2) + pow(z,2),1.5))*\
           (-(vx*z) + x*vz) + \
          (-(vy*z) + y*vz)*\
           (-(pow(vy,2)/mu) + \
             pow(y,2)/\
              pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)\
              - 1/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ) + vz*(-((vy*(x*vx + y*vy + z*vz))/\
                mu) + y*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)))/\
        (np.sqrt(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2))) - \
       (vz*(-(vy*z) + y*vz)*\
          ((-(vx*z) + x*vz)*\
             (-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu)) + \
            (-(vy*z) + y*vz)*\
             (-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu))))/\
        (pow(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2),1.5)*\
          np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2))))/\
     np.sqrt(1 - pow((-(vx*z) + x*vz)*\
           (-((vx*(x*vx + y*vy + z*vz))/mu) + \
             x*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)) + \
          (-(vy*z) + y*vz)*\
           (-((vy*(x*vx + y*vy + z*vz))/mu) +\
             y*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)),2)/\
        ((pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          (pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2)))));

    jacobian[4, 2] = \
    -((-(((-(vx*z) + x*vz)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             (-(vy*z) + y*vz)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)))*\
           (2*((x*z)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - (vx*vz)/mu)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             2*((y*z)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - (vy*vz)/mu)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             2*(pow(z,2)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - 1/\
                 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                pow(vz,2)/mu + \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)*\
              (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu))))/\
        (2.*np.sqrt(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          pow(pow(-((vx*(x*vx + y*vy + z*vz))/\
                 mu) + x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2),1.5)) + \
       ((-(vx*z) + x*vz)*\
           ((x*z)/\
              pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)\
              - (vx*vz)/mu) + \
          (-(vy*z) + y*vz)*\
           ((y*z)/\
              pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)\
              - (vy*vz)/mu) - \
          vx*(-((vx*(x*vx + y*vy + z*vz))/mu) + \
             x*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)) - \
          vy*(-((vy*(x*vx + y*vy + z*vz))/mu) + \
             y*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)))/\
        (np.sqrt(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2))) - \
       ((2*vx*(vx*z - x*vz) - \
            2*vy*(-(vy*z) + y*vz))*\
          ((-(vx*z) + x*vz)*\
             (-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu)) + \
            (-(vy*z) + y*vz)*\
             (-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu))))/\
        (2.*pow(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2),1.5)*\
          np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2))))/\
     np.sqrt(1 - pow((-(vx*z) + x*vz)*\
           (-((vx*(x*vx + y*vy + z*vz))/mu) + \
             x*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)) + \
          (-(vy*z) + y*vz)*\
           (-((vy*(x*vx + y*vy + z*vz))/mu) + \
             y*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)),2)/\
        ((pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          (pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2)))));

    jacobian[4, 3] =                                                         \
    -((-(((-(vx*z) + x*vz)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             (-(vy*z) + y*vz)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)))*\
           (2*((x*vx)/mu - (x*vx + y*vy + z*vz)/mu)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             2*((2*vx*y)/mu - (x*vy)/mu)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             2*((2*vx*z)/mu - (x*vz)/mu)*\
              (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu))))/\
        (2.*np.sqrt(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          pow(pow(-((vx*(x*vx + y*vy + z*vz))/\
                 mu) + x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2),1.5)) + \
       (((2*vx*y)/mu - (x*vy)/mu)*(-(vy*z) + y*vz) + \
          (-(vx*z) + x*vz)*\
           ((x*vx)/mu - (x*vx + y*vy + z*vz)/mu) - \
          z*(-((vx*(x*vx + y*vy + z*vz))/mu) + \
             x*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)))/\
        (np.sqrt(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2))) - \
       (z*(vx*z - x*vz)*\
          ((-(vx*z) + x*vz)*\
             (-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu)) + \
            (-(vy*z) + y*vz)*\
             (-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu))))/\
        (pow(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2),1.5)*\
          np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2))))/\
     np.sqrt(1 - pow((-(vx*z) + x*vz)*\
           (-((vx*(x*vx + y*vy + z*vz))/mu) + \
             x*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)) + \
          (-(vy*z) + y*vz)*\
           (-((vy*(x*vx + y*vy + z*vz))/mu) + \
             y*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)),2)/\
        ((pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          (pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2)))));

    jacobian[4, 4] =  \
    -((-(((-(vx*z) + x*vz)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             (-(vy*z) + y*vz)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)))*\
           (2*(-((vx*y)/mu) + (2*x*vy)/mu)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             2*((y*vy)/mu - (x*vx + y*vy + z*vz)/mu)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             2*((2*vy*z)/mu - (y*vz)/mu)*\
              (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu))))/\
        (2.*np.sqrt(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          pow(pow(-((vx*(x*vx + y*vy + z*vz))/\
                 mu) + x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2),1.5)) + \
       ((-((vx*y)/mu) + (2*x*vy)/mu)*\
           (-(vx*z) + x*vz) + \
          (-(vy*z) + y*vz)*\
           ((y*vy)/mu - (x*vx + y*vy + z*vz)/mu) - \
          z*(-((vy*(x*vx + y*vy + z*vz))/mu) + \
             y*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)))/\
        (np.sqrt(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2))) + \
       (z*(-(vy*z) + y*vz)*\
          ((-(vx*z) + x*vz)*\
             (-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu)) + \
            (-(vy*z) + y*vz)*\
             (-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu))))/\
        (pow(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2),1.5)*\
          np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2))))/\
     np.sqrt(1 - pow((-(vx*z) + x*vz)*\
           (-((vx*(x*vx + y*vy + z*vz))/mu) + \
             x*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)) + \
          (-(vy*z) + y*vz)*\
           (-((vy*(x*vx + y*vy + z*vz))/mu) + \
             y*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)),2)/\
        ((pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          (pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2)))));

    jacobian[4, 5] = \
    -((-(((-(vx*z) + x*vz)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             (-(vy*z) + y*vz)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)))*\
           (2*(-((vx*z)/mu) + (2*x*vz)/mu)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             2*(-((vy*z)/mu) + (2*y*vz)/mu)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) + \
             2*((z*vz)/mu - (x*vx + y*vy + z*vz)/mu)*\
              (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu))))/\
        (2.*np.sqrt(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          pow(pow(-((vx*(x*vx + y*vy + z*vz))/\
                 mu) + x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2),1.5)) + \
       ((-(vx*z) + x*vz)*\
           (-((vx*z)/mu) + (2*x*vz)/mu) + \
          (-(vy*z) + y*vz)*\
           (-((vy*z)/mu) + (2*y*vz)/mu) + \
          x*(-((vx*(x*vx + y*vy + z*vz))/mu) + \
             x*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)) + \
          y*(-((vy*(x*vx + y*vy + z*vz))/mu) + \
             y*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)))/\
        (np.sqrt(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2))) - \
       ((-2*x*(vx*z - x*vz) + 2*y*(-(vy*z) + y*vz))*\
          ((-(vx*z) + x*vz)*\
             (-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu)) + \
            (-(vy*z) + y*vz)*\
             (-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu))))/\
        (2.*pow(pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2),1.5)*\
          np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2))))/\
     np.sqrt(1 - pow((-(vx*z) + x*vz)*\
           (-((vx*(x*vx + y*vy + z*vz))/mu) + \
             x*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)) + \
          (-(vy*z) + y*vz)*\
           (-((vy*(x*vx + y*vy + z*vz))/mu) + \
             y*(-(1/\
                   np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                 + (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)),2)/\
        ((pow(vx*z - x*vz,2) + \
            pow(-(vy*z) + y*vz,2))*\
          (pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
              x*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
              y*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2) + \
            pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
              z*(-(1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))\
                   + (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu),2)))));


    # partials of T (time from periapse) now */

    jacobian[5, 0] = \
    -((((x*vx + y*vy + z*vz)*\
           np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             )*(-2*(-(pow(vx,2)/mu) + \
                pow(x,2)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - 1/\
                 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) - \
             2*(-((vx*vy)/mu) + \
                (x*y)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5))*(-((vy*(x*vx + y*vy + z*vz))/\
                   mu) + y*\
                 (-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) - \
             2*((x*z)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - (vx*vz)/mu)*\
              (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu))))/\
         (2.*mu*pow(1 - \
             pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2),1.5)) + \
        (2*x*(x*vx + y*vy + z*vz)*\
           np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2)))/\
         (mu*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*\
           np.sqrt(1 - pow(-((vx*(x*vx + y*vy + z*vz))/\
                  mu) + x*\
                (-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                  (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) - \
        ((x*vx + y*vy + z*vz)*\
           (2*vy*(-(vx*y) + x*vy) - \
             2*vz*(vx*z - x*vz))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ))/\
         (2.*mu*np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           np.sqrt(1 - pow(-((vx*(x*vx + y*vy + z*vz))/\
                  mu) + x*\
                (-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                  (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) - \
        (vx*np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ))/\
         (mu*np.sqrt(1 - pow(-((vx*\
                    (x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) + \
        (-((x*vx + y*vy + z*vz)*\
               np.sqrt(pow(-(vx*y) + x*vy,2) + \
                 pow(vx*z - x*vz,2) + \
                 pow(-(vy*z) + y*vz,2))*\
               (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                 (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu)*\
               (2*(-(pow(vx,2)/mu) + \
                    pow(x,2)/\
                     pow(pow(x,2) + pow(y,2) + \
                      pow(z,2),1.5) - \
                    1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                      + (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)*\
                  (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                    x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) + \
                 2*(-((vx*vy)/mu) + \
                    (x*y)/\
                     pow(pow(x,2) + pow(y,2) + \
                      pow(z,2),1.5))*\
                  (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                    y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) + \
                 2*((x*z)/\
                     pow(pow(x,2) + pow(y,2) + \
                      pow(z,2),1.5) - (vx*vz)/mu)*\
                  (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                    z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu))))/\
            (2.*mu*np.sqrt(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              pow(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2),1.5)) - \
           ((x*vx + y*vy + z*vz)*\
              np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)*\
              (-2*(-(pow(vx,2)/mu) + \
                   pow(x,2)/\
                    pow(pow(x,2) + pow(y,2) + pow(z,2),\
                     1.5) - \
                   1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                    + (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)*\
                 (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                   x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) - \
                2*(-((vx*vy)/mu) + \
                   (x*y)/\
                    pow(pow(x,2) + pow(y,2) + pow(z,2),\
                     1.5))*\
                 (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                   y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) - \
                2*((x*z)/\
                    pow(pow(x,2) + pow(y,2) + pow(z,2),\
                     1.5) - (vx*vz)/mu)*\
                 (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                   z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu))))/\
            (2.*mu*pow(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2),1.5)*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) - \
           (2*x*(x*vx + y*vy + z*vz)*\
              np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2)))/\
            (mu*pow(pow(x,2) + pow(y,2) + pow(z,2),\
               1.5)*np.sqrt(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) + \
           ((x*vx + y*vy + z*vz)*\
              (2*vy*(-(vx*y) + x*vy) - \
                2*vz*(vx*z - x*vz))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu))/\
            (2.*mu*np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              np.sqrt(1 - pow(-((vx*\
                      (x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) + \
           (vx*np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu))/\
            (mu*np.sqrt(1 - pow(-((vx*\
                       (x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))))/\
         np.sqrt(1 - (pow(x*vx + y*vy + z*vz,2)*\
              (pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              pow(2/\
                 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu,2))/\
            (pow(mu,2)*(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              (pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2)))))/\
      np.sqrt(mu*pow(2/\
           np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
          (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,3)\
        )) - (3*mu*x*pow(2/\
         np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
        (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,2)*\
      (-(((x*vx + y*vy + z*vz)*\
             np.sqrt(pow(-(vx*y) + x*vy,2) + \
               pow(vx*z - x*vz,2) + \
               pow(-(vy*z) + y*vz,2))*\
             (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
               (pow(vx,2) + pow(vy,2) + pow(vz,2))/\
                mu))/\
           (mu*np.sqrt(1 - pow(-((vx*\
                      (x*vx + y*vy + z*vz))/mu) + \
                 x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
               pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                 y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
               pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                 z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2)))) + \
        np.arcsin(((x*vx + y*vy + z*vz)*\
            np.sqrt(pow(-(vx*y) + x*vy,2) + \
              pow(vx*z - x*vz,2) + \
              pow(-(vy*z) + y*vz,2))*\
            (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
              (pow(vx,2) + pow(vy,2) + pow(vz,2))/\
               mu))/\
          (mu*np.sqrt(1 - pow(-((vx*\
                     (x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) - \
              pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) - \
              pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2))*\
            np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                   mu) + x*\
                 (-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) + \
              pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) + \
              pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2))))))/\
    (pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*\
      pow(mu*pow(2/\
           np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
          (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,3)\
        ,1.5));

    jacobian[5, 1] = \
    -((((x*vx + y*vy + z*vz)*\
           np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             )*(-2*(-((vx*vy)/mu) + \
                (x*y)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5))*(-((vx*(x*vx + y*vy + z*vz))/\
                   mu) + x*\
                 (-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) - \
             2*(-(pow(vy,2)/mu) + \
                pow(y,2)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - 1/\
                 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) - \
             2*((y*z)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - (vy*vz)/mu)*\
              (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu))))/\
         (2.*mu*pow(1 - \
             pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2),1.5)) + \
        (2*y*(x*vx + y*vy + z*vz)*\
           np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2)))/\
         (mu*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*\
           np.sqrt(1 - pow(-((vx*(x*vx + y*vy + z*vz))/\
                  mu) + x*\
                (-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                  (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) - \
        ((x*vx + y*vy + z*vz)*\
           (-2*vx*(-(vx*y) + x*vy) + \
             2*vz*(-(vy*z) + y*vz))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ))/\
         (2.*mu*np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           np.sqrt(1 - pow(-((vx*(x*vx + y*vy + z*vz))/\
                  mu) + x*\
                (-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                  (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) - \
        (vy*np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ))/\
         (mu*np.sqrt(1 - pow(-((vx*\
                    (x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) + \
        (-((x*vx + y*vy + z*vz)*\
               np.sqrt(pow(-(vx*y) + x*vy,2) + \
                 pow(vx*z - x*vz,2) + \
                 pow(-(vy*z) + y*vz,2))*\
               (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                 (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu)*\
               (2*(-((vx*vy)/mu) + \
                    (x*y)/\
                     pow(pow(x,2) + pow(y,2) + \
                      pow(z,2),1.5))*\
                  (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                    x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) + \
                 2*(-(pow(vy,2)/mu) + \
                    pow(y,2)/\
                     pow(pow(x,2) + pow(y,2) + \
                      pow(z,2),1.5) - \
                    1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                      + (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)*\
                  (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                    y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) + \
                 2*((y*z)/\
                     pow(pow(x,2) + pow(y,2) + \
                      pow(z,2),1.5) - (vy*vz)/mu)*\
                  (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                    z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu))))/\
            (2.*mu*np.sqrt(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              pow(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2),1.5)) - \
           ((x*vx + y*vy + z*vz)*\
              np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)*\
              (-2*(-((vx*vy)/mu) + \
                   (x*y)/\
                    pow(pow(x,2) + pow(y,2) + pow(z,2),\
                     1.5))*\
                 (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                   x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) - \
                2*(-(pow(vy,2)/mu) + \
                   pow(y,2)/\
                    pow(pow(x,2) + pow(y,2) + pow(z,2),\
                     1.5) - \
                   1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                    + (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)*\
                 (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                   y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) - \
                2*((y*z)/\
                    pow(pow(x,2) + pow(y,2) + pow(z,2),\
                     1.5) - (vy*vz)/mu)*\
                 (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                   z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu))))/\
            (2.*mu*pow(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2),1.5)*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) - \
           (2*y*(x*vx + y*vy + z*vz)*\
              np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2)))/\
            (mu*pow(pow(x,2) + pow(y,2) + pow(z,2),\
               1.5)*np.sqrt(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) + \
           ((x*vx + y*vy + z*vz)*\
              (-2*vx*(-(vx*y) + x*vy) + \
                2*vz*(-(vy*z) + y*vz))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu))/\
            (2.*mu*np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              np.sqrt(1 - pow(-((vx*\
                      (x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) + \
           (vy*np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu))/\
            (mu*np.sqrt(1 - pow(-((vx*\
                       (x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))))/\
         np.sqrt(1 - (pow(x*vx + y*vy + z*vz,2)*\
              (pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              pow(2/\
                 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu,2))/\
            (pow(mu,2)*(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              (pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2)))))/\
      np.sqrt(mu*pow(2/\
           np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
          (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,3)\
        )) - (3*mu*y*pow(2/\
         np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
        (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,2)*\
      (-(((x*vx + y*vy + z*vz)*\
             np.sqrt(pow(-(vx*y) + x*vy,2) + \
               pow(vx*z - x*vz,2) + \
               pow(-(vy*z) + y*vz,2))*\
             (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
               (pow(vx,2) + pow(vy,2) + pow(vz,2))/\
                mu))/\
           (mu*np.sqrt(1 - pow(-((vx*\
                      (x*vx + y*vy + z*vz))/mu) + \
                 x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
               pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                 y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
               pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                 z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2)))) + \
        np.arcsin(((x*vx + y*vy + z*vz)*\
            np.sqrt(pow(-(vx*y) + x*vy,2) + \
              pow(vx*z - x*vz,2) + \
              pow(-(vy*z) + y*vz,2))*\
            (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
              (pow(vx,2) + pow(vy,2) + pow(vz,2))/\
               mu))/\
          (mu*np.sqrt(1 - pow(-((vx*\
                     (x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) - \
              pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) - \
              pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2))*\
            np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                   mu) + x*\
                 (-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) + \
              pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) + \
              pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2))))))/\
    (pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*\
      pow(mu*pow(2/\
           np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
          (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,3)\
        ,1.5));

    jacobian[5, 2] =                                             \
    -((((x*vx + y*vy + z*vz)*\
           np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             )*(-2*((x*z)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - (vx*vz)/mu)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) - \
             2*((y*z)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - (vy*vz)/mu)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) - \
             2*(pow(z,2)/\
                 pow(pow(x,2) + pow(y,2) + pow(z,2),\
                  1.5) - 1/\
                 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                pow(vz,2)/mu + \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)*\
              (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu))))/\
         (2.*mu*pow(1 - \
             pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2),1.5)) + \
        (2*z*(x*vx + y*vy + z*vz)*\
           np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2)))/\
         (mu*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*\
           np.sqrt(1 - pow(-((vx*(x*vx + y*vy + z*vz))/\
                  mu) + x*\
                (-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                  (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) - \
        ((x*vx + y*vy + z*vz)*\
           (2*vx*(vx*z - x*vz) - \
             2*vy*(-(vy*z) + y*vz))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ))/\
         (2.*mu*np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           np.sqrt(1 - pow(-((vx*(x*vx + y*vy + z*vz))/\
                  mu) + x*\
                (-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                  (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) - \
        (vz*np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ))/\
         (mu*np.sqrt(1 - pow(-((vx*\
                    (x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) + \
        (-((x*vx + y*vy + z*vz)*\
               np.sqrt(pow(-(vx*y) + x*vy,2) + \
                 pow(vx*z - x*vz,2) + \
                 pow(-(vy*z) + y*vz,2))*\
               (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                 (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu)*\
               (2*((x*z)/\
                     pow(pow(x,2) + pow(y,2) + \
                      pow(z,2),1.5) - (vx*vz)/mu)*\
                  (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                    x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) + \
                 2*((y*z)/\
                     pow(pow(x,2) + pow(y,2) + \
                      pow(z,2),1.5) - (vy*vz)/mu)*\
                  (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                    y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) + \
                 2*(pow(z,2)/\
                     pow(pow(x,2) + pow(y,2) + \
                      pow(z,2),1.5) - \
                    1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                      - pow(vz,2)/mu + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)*\
                  (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                    z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu))))/\
            (2.*mu*np.sqrt(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              pow(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2),1.5)) - \
           ((x*vx + y*vy + z*vz)*\
              np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)*\
              (-2*((x*z)/\
                    pow(pow(x,2) + pow(y,2) + pow(z,2),\
                     1.5) - (vx*vz)/mu)*\
                 (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                   x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) - \
                2*((y*z)/\
                    pow(pow(x,2) + pow(y,2) + pow(z,2),\
                     1.5) - (vy*vz)/mu)*\
                 (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                   y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) - \
                2*(pow(z,2)/\
                    pow(pow(x,2) + pow(y,2) + pow(z,2),\
                     1.5) - \
                   1/\
                    np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                    - pow(vz,2)/mu + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)*\
                 (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                   z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu))))/\
            (2.*mu*pow(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2),1.5)*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) - \
           (2*z*(x*vx + y*vy + z*vz)*\
              np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2)))/\
            (mu*pow(pow(x,2) + pow(y,2) + pow(z,2),\
               1.5)*np.sqrt(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) + \
           ((x*vx + y*vy + z*vz)*\
              (2*vx*(vx*z - x*vz) - \
                2*vy*(-(vy*z) + y*vz))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu))/\
            (2.*mu*np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              np.sqrt(1 - pow(-((vx*\
                      (x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) + \
           (vz*np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu))/\
            (mu*np.sqrt(1 - pow(-((vx*\
                       (x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))))/\
         np.sqrt(1 - (pow(x*vx + y*vy + z*vz,2)*\
              (pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              pow(2/\
                 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu,2))/\
            (pow(mu,2)*(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              (pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2)))))/\
      np.sqrt(mu*pow(2/\
           np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
          (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,3)\
        )) - (3*mu*z*pow(2/\
         np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
        (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,2)*\
      (-(((x*vx + y*vy + z*vz)*\
             np.sqrt(pow(-(vx*y) + x*vy,2) + \
               pow(vx*z - x*vz,2) + \
               pow(-(vy*z) + y*vz,2))*\
             (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
               (pow(vx,2) + pow(vy,2) + pow(vz,2))/\
                mu))/\
           (mu*np.sqrt(1 - pow(-((vx*\
                      (x*vx + y*vy + z*vz))/mu) + \
                 x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
               pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                 y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
               pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                 z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2)))) + \
        np.arcsin(((x*vx + y*vy + z*vz)*\
            np.sqrt(pow(-(vx*y) + x*vy,2) + \
              pow(vx*z - x*vz,2) + \
              pow(-(vy*z) + y*vz,2))*\
            (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
              (pow(vx,2) + pow(vy,2) + pow(vz,2))/\
               mu))/\
          (mu*np.sqrt(1 - pow(-((vx*\
                     (x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) - \
              pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) - \
              pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2))*\
            np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                   mu) + x*\
                 (-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) + \
              pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) + \
              pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2))))))/\
    (pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*\
      pow(mu*pow(2/\
           np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
          (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,3)\
        ,1.5));

    jacobian[5, 3] = \
    -((((x*vx + y*vy + z*vz)*\
           np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             )*(-2*((x*vx)/mu - \
                (x*vx + y*vy + z*vz)/mu)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) - \
             2*((2*vx*y)/mu - (x*vy)/mu)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) - \
             2*((2*vx*z)/mu - (x*vz)/mu)*\
              (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu))))/\
         (2.*mu*pow(1 - \
             pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2),1.5)) + \
        (2*vx*(x*vx + y*vy + z*vz)*\
           np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2)))/\
         (pow(mu,2)*np.sqrt(1 - \
             pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) - \
        ((x*vx + y*vy + z*vz)*\
           (-2*y*(-(vx*y) + x*vy) + 2*z*(vx*z - x*vz))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ))/\
         (2.*mu*np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           np.sqrt(1 - pow(-((vx*(x*vx + y*vy + z*vz))/\
                  mu) + x*\
                (-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                  (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) - \
        (x*np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ))/\
         (mu*np.sqrt(1 - pow(-((vx*\
                    (x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) + \
        (-((x*vx + y*vy + z*vz)*\
               np.sqrt(pow(-(vx*y) + x*vy,2) + \
                 pow(vx*z - x*vz,2) + \
                 pow(-(vy*z) + y*vz,2))*\
               (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                 (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu)*\
               (2*((x*vx)/mu - \
                    (x*vx + y*vy + z*vz)/mu)*\
                  (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                    x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) + \
                 2*((2*vx*y)/mu - (x*vy)/mu)*\
                  (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                    y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) + \
                 2*((2*vx*z)/mu - (x*vz)/mu)*\
                  (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                    z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu))))/\
            (2.*mu*np.sqrt(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              pow(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2),1.5)) - \
           ((x*vx + y*vy + z*vz)*\
              np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)*\
              (-2*((x*vx)/mu - \
                   (x*vx + y*vy + z*vz)/mu)*\
                 (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                   x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) - \
                2*((2*vx*y)/mu - (x*vy)/mu)*\
                 (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                   y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) - \
                2*((2*vx*z)/mu - (x*vz)/mu)*\
                 (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                   z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu))))/\
            (2.*mu*pow(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2),1.5)*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) - \
           (2*vx*(x*vx + y*vy + z*vz)*\
              np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2)))/\
            (pow(mu,2)*np.sqrt(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) + \
           ((x*vx + y*vy + z*vz)*\
              (-2*y*(-(vx*y) + x*vy) + \
                2*z*(vx*z - x*vz))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu))/\
            (2.*mu*np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              np.sqrt(1 - pow(-((vx*\
                      (x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) + \
           (x*np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu))/\
            (mu*np.sqrt(1 - pow(-((vx*\
                       (x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))))/\
         np.sqrt(1 - (pow(x*vx + y*vy + z*vz,2)*\
              (pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              pow(2/\
                 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu,2))/\
            (pow(mu,2)*(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              (pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2)))))/\
      np.sqrt(mu*pow(2/\
           np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
          (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,3)\
        )) - (3*vx*pow(2/\
         np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
        (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,2)*\
      (-(((x*vx + y*vy + z*vz)*\
             np.sqrt(pow(-(vx*y) + x*vy,2) + \
               pow(vx*z - x*vz,2) + \
               pow(-(vy*z) + y*vz,2))*\
             (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
               (pow(vx,2) + pow(vy,2) + pow(vz,2))/\
                mu))/\
           (mu*np.sqrt(1 - pow(-((vx*\
                      (x*vx + y*vy + z*vz))/mu) + \
                 x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
               pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                 y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
               pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                 z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2)))) + \
        np.arcsin(((x*vx + y*vy + z*vz)*\
            np.sqrt(pow(-(vx*y) + x*vy,2) + \
              pow(vx*z - x*vz,2) + \
              pow(-(vy*z) + y*vz,2))*\
            (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
              (pow(vx,2) + pow(vy,2) + pow(vz,2))/\
               mu))/\
          (mu*np.sqrt(1 - pow(-((vx*\
                     (x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) - \
              pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) - \
              pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2))*\
            np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                   mu) + x*\
                 (-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) + \
              pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) + \
              pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2))))))/\
    pow(mu*pow(2/\
         np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
        (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,3),\
     1.5);


    jacobian[5, 4] =  \
    -((((x*vx + y*vy + z*vz)*\
           np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             )*(-2*(-((vx*y)/mu) + (2*x*vy)/mu)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) - \
             2*((y*vy)/mu - (x*vx + y*vy + z*vz)/mu)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) - \
             2*((2*vy*z)/mu - (y*vz)/mu)*\
              (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu))))/\
         (2.*mu*pow(1 - \
             pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2),1.5)) + \
        (2*vy*(x*vx + y*vy + z*vz)*\
           np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2)))/\
         (pow(mu,2)*np.sqrt(1 - \
             pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) - \
        ((x*vx + y*vy + z*vz)*\
           (2*x*(-(vx*y) + x*vy) - \
             2*z*(-(vy*z) + y*vz))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ))/\
         (2.*mu*np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           np.sqrt(1 - pow(-((vx*(x*vx + y*vy + z*vz))/\
                  mu) + x*\
                (-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                  (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) - \
        (y*np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ))/\
         (mu*np.sqrt(1 - pow(-((vx*\
                    (x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) + \
        (-((x*vx + y*vy + z*vz)*\
               np.sqrt(pow(-(vx*y) + x*vy,2) + \
                 pow(vx*z - x*vz,2) + \
                 pow(-(vy*z) + y*vz,2))*\
               (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                 (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu)*\
               (2*(-((vx*y)/mu) + (2*x*vy)/mu)*\
                  (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                    x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) + \
                 2*((y*vy)/mu - \
                    (x*vx + y*vy + z*vz)/mu)*\
                  (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                    y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) + \
                 2*((2*vy*z)/mu - (y*vz)/mu)*\
                  (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                    z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu))))/\
            (2.*mu*np.sqrt(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              pow(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2),1.5)) - \
           ((x*vx + y*vy + z*vz)*\
              np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)*\
              (-2*(-((vx*y)/mu) + (2*x*vy)/mu)*\
                 (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                   x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) - \
                2*((y*vy)/mu - \
                   (x*vx + y*vy + z*vz)/mu)*\
                 (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                   y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) - \
                2*((2*vy*z)/mu - (y*vz)/mu)*\
                 (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                   z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu))))/\
            (2.*mu*pow(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2),1.5)*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) - \
           (2*vy*(x*vx + y*vy + z*vz)*\
              np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2)))/\
            (pow(mu,2)*np.sqrt(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) + \
           ((x*vx + y*vy + z*vz)*\
              (2*x*(-(vx*y) + x*vy) - \
                2*z*(-(vy*z) + y*vz))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu))/\
            (2.*mu*np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              np.sqrt(1 - pow(-((vx*\
                      (x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) + \
           (y*np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu))/\
            (mu*np.sqrt(1 - pow(-((vx*\
                       (x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))))/\
         np.sqrt(1 - (pow(x*vx + y*vy + z*vz,2)*\
              (pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              pow(2/\
                 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu,2))/\
            (pow(mu,2)*(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              (pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2)))))/\
      np.sqrt(mu*pow(2/\
           np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
          (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,3)\
        )) - (3*vy*pow(2/\
         np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
        (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,2)*\
      (-(((x*vx + y*vy + z*vz)*\
             np.sqrt(pow(-(vx*y) + x*vy,2) + \
               pow(vx*z - x*vz,2) + \
               pow(-(vy*z) + y*vz,2))*\
             (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
               (pow(vx,2) + pow(vy,2) + pow(vz,2))/\
                mu))/\
           (mu*np.sqrt(1 - pow(-((vx*\
                      (x*vx + y*vy + z*vz))/mu) + \
                 x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
               pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                 y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
               pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                 z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2)))) + \
        np.arcsin(((x*vx + y*vy + z*vz)*\
            np.sqrt(pow(-(vx*y) + x*vy,2) + \
              pow(vx*z - x*vz,2) + \
              pow(-(vy*z) + y*vz,2))*\
            (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
              (pow(vx,2) + pow(vy,2) + pow(vz,2))/\
               mu))/\
          (mu*np.sqrt(1 - pow(-((vx*\
                     (x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) - \
              pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) - \
              pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2))*\
            np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                   mu) + x*\
                 (-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) + \
              pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) + \
              pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2))))))/\
    pow(mu*pow(2/\
         np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
        (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,3),\
     1.5);

    jacobian[5, 5] =                                              \
    -((((x*vx + y*vy + z*vz)*\
           np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             )*(-2*(-((vx*z)/mu) + (2*x*vz)/mu)*\
              (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) - \
             2*(-((vy*z)/mu) + (2*y*vz)/mu)*\
              (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu)) - \
             2*((z*vz)/mu - (x*vx + y*vy + z*vz)/mu)*\
              (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu))))/\
         (2.*mu*pow(1 - \
             pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2),1.5)) + \
        (2*vz*(x*vx + y*vy + z*vz)*\
           np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2)))/\
         (pow(mu,2)*np.sqrt(1 - \
             pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) - \
        ((x*vx + y*vy + z*vz)*\
           (-2*x*(vx*z - x*vz) + 2*y*(-(vy*z) + y*vz))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ))/\
         (2.*mu*np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           np.sqrt(1 - pow(-((vx*(x*vx + y*vy + z*vz))/\
                  mu) + x*\
                (-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                  (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) - \
        (z*np.sqrt(pow(-(vx*y) + x*vy,2) + \
             pow(vx*z - x*vz,2) + \
             pow(-(vy*z) + y*vz,2))*\
           (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
             (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu\
             ))/\
         (mu*np.sqrt(1 - pow(-((vx*\
                    (x*vx + y*vy + z*vz))/mu) + \
               x*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
               y*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2) - \
             pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
               z*(-(1/\
                     np.sqrt(pow(x,2) + pow(y,2) + pow(z,2))\
                     ) + (pow(vx,2) + pow(vy,2) + \
                     pow(vz,2))/mu),2))) + \
        (-((x*vx + y*vy + z*vz)*\
               np.sqrt(pow(-(vx*y) + x*vy,2) + \
                 pow(vx*z - x*vz,2) + \
                 pow(-(vy*z) + y*vz,2))*\
               (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                 (pow(vx,2) + pow(vy,2) + \
                    pow(vz,2))/mu)*\
               (2*(-((vx*z)/mu) + (2*x*vz)/mu)*\
                  (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                    x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) + \
                 2*(-((vy*z)/mu) + (2*y*vz)/mu)*\
                  (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                    y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) + \
                 2*((z*vz)/mu - \
                    (x*vx + y*vy + z*vz)/mu)*\
                  (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                    z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                       (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu))))/\
            (2.*mu*np.sqrt(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              pow(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2),1.5)) - \
           ((x*vx + y*vy + z*vz)*\
              np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu)*\
              (-2*(-((vx*z)/mu) + (2*x*vz)/mu)*\
                 (-((vx*(x*vx + y*vy + z*vz))/mu) + \
                   x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) - \
                2*(-((vy*z)/mu) + (2*y*vz)/mu)*\
                 (-((vy*(x*vx + y*vy + z*vz))/mu) + \
                   y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu)) - \
                2*((z*vz)/mu - \
                   (x*vx + y*vy + z*vz)/mu)*\
                 (-((vz*(x*vx + y*vy + z*vz))/mu) + \
                   z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                      (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu))))/\
            (2.*mu*pow(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2),1.5)*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) - \
           (2*vz*(x*vx + y*vy + z*vz)*\
              np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2)))/\
            (pow(mu,2)*np.sqrt(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) + \
           ((x*vx + y*vy + z*vz)*\
              (-2*x*(vx*z - x*vz) + \
                2*y*(-(vy*z) + y*vz))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu))/\
            (2.*mu*np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              np.sqrt(1 - pow(-((vx*\
                      (x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))) + \
           (z*np.sqrt(pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu))/\
            (mu*np.sqrt(1 - pow(-((vx*\
                       (x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))))/\
         np.sqrt(1 - (pow(x*vx + y*vy + z*vz,2)*\
              (pow(-(vx*y) + x*vy,2) + \
                pow(vx*z - x*vz,2) + \
                pow(-(vy*z) + y*vz,2))*\
              pow(2/\
                 np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
                (pow(vx,2) + pow(vy,2) + \
                   pow(vz,2))/mu,2))/\
            (pow(mu,2)*(1 - \
                pow(-((vx*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2))*\
              (pow(-((vx*(x*vx + y*vy + z*vz))/mu) + \
                  x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vy*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) + \
                pow(-((vz*(x*vx + y*vy + z*vz))/\
                     mu) + \
                  z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                     (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2)))))/\
      np.sqrt(mu*pow(2/\
           np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
          (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,3)\
        )) - (3*vz*pow(2/\
         np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
        (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,2)*\
      (-(((x*vx + y*vy + z*vz)*\
             np.sqrt(pow(-(vx*y) + x*vy,2) + \
               pow(vx*z - x*vz,2) + \
               pow(-(vy*z) + y*vz,2))*\
             (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
               (pow(vx,2) + pow(vy,2) + pow(vz,2))/\
                mu))/\
           (mu*np.sqrt(1 - pow(-((vx*\
                      (x*vx + y*vy + z*vz))/mu) + \
                 x*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
               pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                 y*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2) - \
               pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                 z*(-(1/\
                       np.sqrt(pow(x,2) + pow(y,2) + \
                       pow(z,2))) + \
                    (pow(vx,2) + pow(vy,2) + \
                       pow(vz,2))/mu),2)))) + \
        np.arcsin(((x*vx + y*vy + z*vz)*\
            np.sqrt(pow(-(vx*y) + x*vy,2) + \
              pow(vx*z - x*vz,2) + \
              pow(-(vy*z) + y*vz,2))*\
            (2/np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
              (pow(vx,2) + pow(vy,2) + pow(vz,2))/\
               mu))/\
          (mu*np.sqrt(1 - pow(-((vx*\
                     (x*vx + y*vy + z*vz))/mu) + \
                x*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) - \
              pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) - \
              pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2))*\
            np.sqrt(pow(-((vx*(x*vx + y*vy + z*vz))/\
                   mu) + x*\
                 (-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) + \
              pow(-((vy*(x*vx + y*vy + z*vz))/mu) + \
                y*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2) + \
              pow(-((vz*(x*vx + y*vy + z*vz))/mu) + \
                z*(-(1/\
                      np.sqrt(pow(x,2) + pow(y,2) + \
                      pow(z,2))) + \
                   (pow(vx,2) + pow(vy,2) + \
                      pow(vz,2))/mu),2))))))/\
    pow(mu*pow(2/\
         np.sqrt(pow(x,2) + pow(y,2) + pow(z,2)) - \
        (pow(vx,2) + pow(vy,2) + pow(vz,2))/mu,3),\
     1.5);

    return jacobian
