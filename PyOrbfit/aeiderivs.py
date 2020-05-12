def aeiderivs(self):
  aei_partials = np.self.Zeros([6, 6])
  # aself.X is del a bself.Y del self.X */
  aei_partials[1][1] = \
    (2*self.X)/(pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5)*\
	   pow(2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
	       (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,2));

  aei_partials[1][2] = \
  (2*self.Y)/(pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5)*\
	    pow(2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
		(pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,2));

  aei_partials[1][3] = \
  (2*self.Z)/(pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5)*\
	    pow(2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
		(pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,2));

  aei_partials[1][4] = \
  (2*self.VX)/(μ*pow(2/\
			 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
			 (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,2));

  aei_partials[1][5] = \
  (2*self.VY)/(μ*pow(2/\
			 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
			 (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,2));

  aei_partials[1][6] = \
    (2*self.VZ)/(μ*pow(2/\
		     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
		     (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,2));


# Partials for e now */

  aei_partials[2][1] = \
    (2*(-(pow(self.VX,2)/μ) + \
        pow(self.X,2)/\
         pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5) - \
        1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) + \
        (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)*\
      (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.X*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
       + 2*(-((self.VX*self.VY)/μ) + \
        (self.X*self.Y)/pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5))*\
      (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.Y*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
       + 2*((self.X*self.Z)/\
         pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5) - \
        (self.VX*self.VZ)/μ)*(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
           μ) + self.Z*(-(1/\
              np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
     )/(2.*np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
            μ) + self.X*(-(1/\
               np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2) + pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
         self.Y*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2) + pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
         self.Z*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2)));

  aei_partials[2][2] = \
    (2*(-((self.VX*self.VY)/μ) + (self.X*self.Y)/\
         pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5))*\
      (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.X*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
       + 2*(-(pow(self.VY,2)/μ) + \
        pow(self.Y,2)/\
         pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5) - \
        1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) + \
        (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)*\
      (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.Y*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
       + 2*((self.Y*self.Z)/\
         pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5) - \
        (self.VY*self.VZ)/μ)*(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
           μ) + self.Z*(-(1/\
              np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
     )/(2.*np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
            μ) + self.X*(-(1/\
               np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2) + pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
         self.Y*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2) + pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
         self.Z*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2)));

  aei_partials[2][3] = \
    (2*((self.X*self.Z)/pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5) - \
        (self.VX*self.VZ)/μ)*(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
           μ) + self.X*(-(1/\
              np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
       + 2*((self.Y*self.Z)/\
         pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5) - \
        (self.VY*self.VZ)/μ)*(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
           μ) + self.Y*(-(1/\
              np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
       + 2*(pow(self.Z,2)/\
         pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5) - \
        1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
        pow(self.VZ,2)/μ + \
        (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)*\
      (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.Z*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
     )/(2.*np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
            μ) + self.X*(-(1/\
               np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2) + pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
         self.Y*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2) + pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
         self.Z*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2)));

  aei_partials[2][4] = \
    (2*((self.X*self.VX)/μ - (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
      (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.X*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
       + 2*((2*self.VX*self.Y)/μ - (self.X*self.VY)/μ)*\
      (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.Y*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
       + 2*((2*self.VX*self.Z)/μ - (self.X*self.VZ)/μ)*\
      (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.Z*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
     )/(2.*np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
            μ) + self.X*(-(1/\
               np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2) + pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
         self.Y*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2) + pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
         self.Z*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2)));

  aei_partials[2][5] = \
    (2*(-((self.VX*self.Y)/μ) + (2*self.X*self.VY)/μ)*\
      (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.X*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
       + 2*((self.Y*self.VY)/μ - (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
      (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.Y*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
       + 2*((2*self.VY*self.Z)/μ - (self.Y*self.VZ)/μ)*\
      (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.Z*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
     )/(2.*np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
            μ) + self.X*(-(1/\
               np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2) + pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
         self.Y*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2) + pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
         self.Z*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2)));

  aei_partials[2][6] = \
    (2*(-((self.VX*self.Z)/μ) + (2*self.X*self.VZ)/μ)*\
      (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.X*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
       + 2*(-((self.VY*self.Z)/μ) + (2*self.Y*self.VZ)/μ)*\
      (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.Y*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
       + 2*((self.Z*self.VZ)/μ - (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
      (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
        self.Z*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
           (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ))\
     )/(2.*np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
            μ) + self.X*(-(1/\
               np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) +\
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2) + pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
         self.Y*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2) + pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
         self.Z*(-(1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))) + \
            (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ)\
         ,2)));

# Partials of i now */

  aei_partials[3][1] = \
    -((-((-(self.VX*self.Y) + self.X*self.VY)*\
           (2*self.VY*(-(self.VX*self.Y) + self.X*self.VY) - \
             2*self.VZ*(self.VX*self.Z - self.X*self.VZ)))/\
        (2.*pow(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
            pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)) + \
       self.VY/np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
          pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2)\
          ))/\
     np.sqrt(1 - pow(-(self.VX*self.Y) + self.X*self.VY,2)/\
        (pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
          pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2)\
          )));

  aei_partials[3][2] = \
    -((-((-(self.VX*self.Y) + self.X*self.VY)*\
           (-2*self.VX*(-(self.VX*self.Y) + self.X*self.VY) + \
             2*self.VZ*(-(self.VY*self.Z) + self.Y*self.VZ)))/\
        (2.*pow(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
            pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)) - \
       self.VX/np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
          pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2)\
          ))/\
     np.sqrt(1 - pow(-(self.VX*self.Y) + self.X*self.VY,2)/\
        (pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
          pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2)\
          )));

  aei_partials[3][3] = \
    ((-(self.VX*self.Y) + self.X*self.VY)*(2*self.VX*(self.VX*self.Z - self.X*self.VZ) - \
       2*self.VY*(-(self.VY*self.Z) + self.Y*self.VZ)))/\
   (2.*pow(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
       pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2),\
      1.5)*np.sqrt(1 - pow(-(self.VX*self.Y) + self.X*self.VY,2)/\
        (pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
          pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2)\
          )));

  aei_partials[3][4] = \
    -((-((-(self.VX*self.Y) + self.X*self.VY)*\
           (-2*self.Y*(-(self.VX*self.Y) + self.X*self.VY) + 2*self.Z*(self.VX*self.Z - self.X*self.VZ)))/
        (2.*pow(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
            pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)) - \
       self.Y/np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
          pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2)\
          ))/\
     np.sqrt(1 - pow(-(self.VX*self.Y) + self.X*self.VY,2)/\
        (pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
          pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2)\
          )));

  aei_partials[3][5] = \
    -((-((-(self.VX*self.Y) + self.X*self.VY)*\
           (2*self.X*(-(self.VX*self.Y) + self.X*self.VY) - 2*self.Z*(-(self.VY*self.Z) + self.Y*self.VZ))\
           )/\
        (2.*pow(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
            pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)) + \
       self.X/np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
          pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2)\
          ))/\
     np.sqrt(1 - pow(-(self.VX*self.Y) + self.X*self.VY,2)/\
        (pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
          pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2)\
          )));

  aei_partials[3][6] = \
    ((-(self.VX*self.Y) + self.X*self.VY)*(-2*self.X*(self.VX*self.Z - self.X*self.VZ) + \
       2*self.Y*(-(self.VY*self.Z) + self.Y*self.VZ)))/\
   (2.*pow(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
       pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2),\
      1.5)*np.sqrt(1 - pow(-(self.VX*self.Y) + self.X*self.VY,2)/\
        (pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
          pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2)\
          )));


# Partials of capital Omega (long of asc node) now */

  aei_partials[4][1] = \
    -(((self.VZ*(self.VX*self.Z - self.X*self.VZ)*(-(self.VX*self.Z) + self.X*self.VZ))/\
        pow(pow(self.VX*self.Z - self.X*self.VZ,2) + \
          pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5) + \
       self.VZ/np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
          pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
     np.sqrt(1 - pow(-(self.VX*self.Z) + self.X*self.VZ,2)/\
        (pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2))\
       ));

  aei_partials[4][2] = \
    (self.VZ*(-(self.VX*self.Z) + self.X*self.VZ)*(-(self.VY*self.Z) + self.Y*self.VZ))/\
   (pow(pow(self.VX*self.Z - self.X*self.VZ,2) + \
       pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)*\
     np.sqrt(1 - pow(-(self.VX*self.Z) + self.X*self.VZ,2)/\
        (pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2))\
       ));

  aei_partials[4][3] = \
    -((-((-(self.VX*self.Z) + self.X*self.VZ)*\
           (2*self.VX*(self.VX*self.Z - self.X*self.VZ) - \
             2*self.VY*(-(self.VY*self.Z) + self.Y*self.VZ)))/\
        (2.*pow(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)) - \
       self.VX/np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
          pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
     np.sqrt(1 - pow(-(self.VX*self.Z) + self.X*self.VZ,2)/\
        (pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2))\
       ));

  aei_partials[4][4] = \
    -((-((self.Z*(self.VX*self.Z - self.X*self.VZ)*(-(self.VX*self.Z) + self.X*self.VZ))/\
          pow(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)) - \
       self.Z/np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
          pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
     np.sqrt(1 - pow(-(self.VX*self.Z) + self.X*self.VZ,2)/\
        (pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2))\
       ));

  aei_partials[4][5] = \
    -((self.Z*(-(self.VX*self.Z) + self.X*self.VZ)*(-(self.VY*self.Z) + self.Y*self.VZ))/\
     (pow(pow(self.VX*self.Z - self.X*self.VZ,2) + \
         pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)*\
       np.sqrt(1 - pow(-(self.VX*self.Z) + self.X*self.VZ,2)/\
          (pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))));

  aei_partials[4][6] = \
    -((-((-(self.VX*self.Z) + self.X*self.VZ)*\
           (-2*self.X*(self.VX*self.Z - self.X*self.VZ) + 2*self.Y*(-(self.VY*self.Z) + self.Y*self.VZ)))/\
        (2.*pow(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)) + \
       self.X/np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
          pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
     np.sqrt(1 - pow(-(self.VX*self.Z) + self.X*self.VZ,2)/\
        (pow(self.VX*self.Z - self.X*self.VZ,2) + pow(-(self.VY*self.Z) + self.Y*self.VZ,2))\
       ));


# partials of small omega (arg of per) now */

  aei_partials[5][1] = \
    -((-(((-(self.VX*self.Z) + self.X*self.VZ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             (-(self.VY*self.Z) + self.Y*self.VZ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)))*\
           (2*(-(pow(self.VX,2)/μ) + \
                pow(self.X,2)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - 1/\
                 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) +\
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             2*(-((self.VX*self.VY)/μ) + \
                (self.X*self.Y)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5))*(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                   μ) + self.Y*\
                 (-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             2*((self.X*self.Z)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - (self.VX*self.VZ)/μ)*\
              (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ))))/\
        (2.*np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          pow(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                 μ) + self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2),1.5)) + \
       ((-((self.VX*self.VY)/μ) + \
             (self.X*self.Y)/\
              pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5))*\
           (-(self.VY*self.Z) + self.Y*self.VZ) + \
          (-(self.VX*self.Z) + self.X*self.VZ)*\
           (-(pow(self.VX,2)/μ) + \
             pow(self.X,2)/\
              pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5)\
              - 1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) + \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ) + self.VZ*(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                μ) + self.X*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)))/\
        (np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2))) + \
       (self.VZ*(self.VX*self.Z - self.X*self.VZ)*\
          ((-(self.VX*self.Z) + self.X*self.VZ)*\
             (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ)) + \
            (-(self.VY*self.Z) + self.Y*self.VZ)*\
             (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ))))/\
        (pow(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)*\
          np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2))))/\
     np.sqrt(1 - pow((-(self.VX*self.Z) + self.X*self.VZ)*\
           (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.X*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)) + \
          (-(self.VY*self.Z) + self.Y*self.VZ)*\
           (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.Y*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)),2)/\
        ((pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          (pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2)))));

  aei_partials[5][2] = \
    -((-(((-(self.VX*self.Z) + self.X*self.VZ)*
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) +\
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow\(self.VY,2) + \
                      pow(self.VZ,2))/μ)\) + \
             (-(self.VY*self.Z) + self.Y*self.VZ)*\
              (-((self.VY*(self.X*self.VX + self.Y\*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow\(self.VY,2) + \
                      pow(self.VZ,2))/μ)\))*\
           (2*(-((self.VX*self.VY)/μ) + \
                (self.X*self.Y)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5)\)*(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                   μ) + self.X*\
                 (-(1/
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             2*(-(pow(self.VY,2)/μ) + \
                pow(self.Y,2)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - 1/\
                 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) + \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             2*((self.Y*self.Z)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - (self.VY*self.VZ)/μ)*\
              (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ))))/\
        (2.*np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          pow(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                 μ) + self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\)\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2),1.5)) + \
       ((-((self.VX*self.VY)/μ) + \
             (self.X*self.Y)/\
              pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5))*\
           (-(self.VX*self.Z) + self.X*self.VZ) + \
          (-(self.VY*self.Z) + self.Y*self.VZ)*\
           (-(pow(self.VY,2)/μ) + \
             pow(self.Y,2)/\
              pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5)\
              - 1/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) + \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ) + self.VZ*(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                μ) + self.Y*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)))/\
        (np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2))) - \
       (self.VZ*(-(self.VY*self.Z) + self.Y*self.VZ)*\
          ((-(self.VX*self.Z) + self.X*self.VZ)*\
             (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ)) + \
            (-(self.VY*self.Z) + self.Y*self.VZ)*\
             (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ))))/\
        (pow(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)*\
          np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2))))/\
     np.sqrt(1 - pow((-(self.VX*self.Z) + self.X*self.VZ)*\
           (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.X*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)) + \
          (-(self.VY*self.Z) + self.Y*self.VZ)*\
           (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) +\
             self.Y*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)),2)/\
        ((pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          (pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2)))));

  aei_partials[5][3] = \
    -((-(((-(self.VX*self.Z) + self.X*self.VZ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             (-(self.VY*self.Z) + self.Y*self.VZ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)))*\
           (2*((self.X*self.Z)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - (self.VX*self.VZ)/μ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             2*((self.Y*self.Z)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - (self.VY*self.VZ)/μ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             2*(pow(self.Z,2)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - 1/\
                 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                pow(self.VZ,2)/μ + \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)*\
              (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ))))/\
        (2.*np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          pow(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                 μ) + self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2),1.5)) + \
       ((-(self.VX*self.Z) + self.X*self.VZ)*\
           ((self.X*self.Z)/\
              pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5)\
              - (self.VX*self.VZ)/μ) + \
          (-(self.VY*self.Z) + self.Y*self.VZ)*\
           ((self.Y*self.Z)/\
              pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5)\
              - (self.VY*self.VZ)/μ) - \
          self.VX*(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.X*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)) - \
          self.VY*(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.Y*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)))/\
        (np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2))) - \
       ((2*self.VX*(self.VX*self.Z - self.X*self.VZ) - \
            2*self.VY*(-(self.VY*self.Z) + self.Y*self.VZ))*\
          ((-(self.VX*self.Z) + self.X*self.VZ)*\
             (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ)) + \
            (-(self.VY*self.Z) + self.Y*self.VZ)*\
             (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ))))/\
        (2.*pow(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)*\
          np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2))))/\
     np.sqrt(1 - pow((-(self.VX*self.Z) + self.X*self.VZ)*\
           (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.X*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)) + \
          (-(self.VY*self.Z) + self.Y*self.VZ)*\
           (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.Y*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)),2)/\
        ((pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          (pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2)))));

  aei_partials[5][4] =                                                         \
    -((-(((-(self.VX*self.Z) + self.X*self.VZ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             (-(self.VY*self.Z) + self.Y*self.VZ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)))*\
           (2*((self.X*self.VX)/μ - (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             2*((2*self.VX*self.Y)/μ - (self.X*self.VY)/μ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             2*((2*self.VX*self.Z)/μ - (self.X*self.VZ)/μ)*\
              (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ))))/\
        (2.*np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          pow(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                 μ) + self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2),1.5)) + \
       (((2*self.VX*self.Y)/μ - (self.X*self.VY)/μ)*(-(self.VY*self.Z) + self.Y*self.VZ) + \
          (-(self.VX*self.Z) + self.X*self.VZ)*\
           ((self.X*self.VX)/μ - (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ) - \
          self.Z*(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.X*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)))/\
        (np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2))) - \
       (self.Z*(self.VX*self.Z - self.X*self.VZ)*\
          ((-(self.VX*self.Z) + self.X*self.VZ)*\
             (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ)) + \
            (-(self.VY*self.Z) + self.Y*self.VZ)*\
             (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ))))/\
        (pow(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)*\
          np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2))))/\
     np.sqrt(1 - pow((-(self.VX*self.Z) + self.X*self.VZ)*\
           (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.X*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)) + \
          (-(self.VY*self.Z) + self.Y*self.VZ)*\
           (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.Y*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)),2)/\
        ((pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          (pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2)))));

  aei_partials[5][5] =  \
    -((-(((-(self.VX*self.Z) + self.X*self.VZ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             (-(self.VY*self.Z) + self.Y*self.VZ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)))*\
           (2*(-((self.VX*self.Y)/μ) + (2*self.X*self.VY)/μ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             2*((self.Y*self.VY)/μ - (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             2*((2*self.VY*self.Z)/μ - (self.Y*self.VZ)/μ)*\
              (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ))))/\
        (2.*np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          pow(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                 μ) + self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2),1.5)) + \
       ((-((self.VX*self.Y)/μ) + (2*self.X*self.VY)/μ)*\
           (-(self.VX*self.Z) + self.X*self.VZ) + \
          (-(self.VY*self.Z) + self.Y*self.VZ)*\
           ((self.Y*self.VY)/μ - (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ) - \
          self.Z*(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.Y*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)))/\
        (np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2))) + \
       (self.Z*(-(self.VY*self.Z) + self.Y*self.VZ)*\
          ((-(self.VX*self.Z) + self.X*self.VZ)*\
             (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ)) + \
            (-(self.VY*self.Z) + self.Y*self.VZ)*\
             (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ))))/\
        (pow(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)*\
          np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2))))/\
     np.sqrt(1 - pow((-(self.VX*self.Z) + self.X*self.VZ)*\
           (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.X*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)) + \
          (-(self.VY*self.Z) + self.Y*self.VZ)*\
           (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.Y*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)),2)/\
        ((pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          (pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2)))));

  aei_partials[5][6] = \
    -((-(((-(self.VX*self.Z) + self.X*self.VZ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             (-(self.VY*self.Z) + self.Y*self.VZ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)))*\
           (2*(-((self.VX*self.Z)/μ) + (2*self.X*self.VZ)/μ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             2*(-((self.VY*self.Z)/μ) + (2*self.Y*self.VZ)/μ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) + \
             2*((self.Z*self.VZ)/μ - (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
              (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ))))/\
        (2.*np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          pow(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                 μ) + self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2),1.5)) + \
       ((-(self.VX*self.Z) + self.X*self.VZ)*\
           (-((self.VX*self.Z)/μ) + (2*self.X*self.VZ)/μ) + \
          (-(self.VY*self.Z) + self.Y*self.VZ)*\
           (-((self.VY*self.Z)/μ) + (2*self.Y*self.VZ)/μ) + \
          self.X*(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.X*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)) + \
          self.Y*(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.Y*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)))/\
        (np.sqrt(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2))) - \
       ((-2*self.X*(self.VX*self.Z - self.X*self.VZ) + 2*self.Y*(-(self.VY*self.Z) + self.Y*self.VZ))*\
          ((-(self.VX*self.Z) + self.X*self.VZ)*\
             (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ)) + \
            (-(self.VY*self.Z) + self.Y*self.VZ)*\
             (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ))))/\
        (2.*pow(pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2),1.5)*\
          np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2))))/\
     np.sqrt(1 - pow((-(self.VX*self.Z) + self.X*self.VZ)*\
           (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.X*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)) + \
          (-(self.VY*self.Z) + self.Y*self.VZ)*\
           (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
             self.Y*(-(1/\
                   np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                 + (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)),2)/\
        ((pow(self.VX*self.Z - self.X*self.VZ,2) + \
            pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
          (pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.X*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Y*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2) + \
            pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
              self.Z*(-(1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)))\
                   + (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ),2)))));


# partials of T (time from periapse) now */

  aei_partials[6][1] = \
    -((((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             )*(-2*(-(pow(self.VX,2)/μ) + \
                pow(self.X,2)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - 1/\
                 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) + \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) - \
             2*(-((self.VX*self.VY)/μ) + \
                (self.X*self.Y)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5))*(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                   μ) + self.Y*\
                 (-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) - \
             2*((self.X*self.Z)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - (self.VX*self.VZ)/μ)*\
              (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ))))/\
         (2.*μ*pow(1 - \
             pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2),1.5)) + \
        (2*self.X*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
         (μ*pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5)*\
           np.sqrt(1 - pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                  μ) + self.X*\
                (-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                  (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) - \
        ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           (2*self.VY*(-(self.VX*self.Y) + self.X*self.VY) - \
             2*self.VZ*(self.VX*self.Z - self.X*self.VZ))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ))/\
         (2.*μ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           np.sqrt(1 - pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                  μ) + self.X*\
                (-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                  (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) - \
        (self.VX*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ))/\
         (μ*np.sqrt(1 - pow(-((self.VX*\
                    (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) + \
        (-((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
               np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                 pow(self.VX*self.Z - self.X*self.VZ,2) + \
                 pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
               (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                 (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ)*\
               (2*(-(pow(self.VX,2)/μ) + \
                    pow(self.X,2)/\
                     pow(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2),1.5) - \
                    1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                      + (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)*\
                  (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) + \
                 2*(-((self.VX*self.VY)/μ) + \
                    (self.X*self.Y)/\
                     pow(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2),1.5))*\
                  (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) + \
                 2*((self.X*self.Z)/\
                     pow(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2),1.5) - (self.VX*self.VZ)/μ)*\
                  (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ))))/\
            (2.*μ*np.sqrt(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              pow(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2),1.5)) - \
           ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)*\
              (-2*(-(pow(self.VX,2)/μ) + \
                   pow(self.X,2)/\
                    pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                     1.5) - \
                   1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                    + (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)*\
                 (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) - \
                2*(-((self.VX*self.VY)/μ) + \
                   (self.X*self.Y)/\
                    pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                     1.5))*\
                 (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) - \
                2*((self.X*self.Z)/\
                    pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                     1.5) - (self.VX*self.VZ)/μ)*\
                 (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ))))/\
            (2.*μ*pow(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2),1.5)*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) - \
           (2*self.X*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
            (μ*pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
               1.5)*np.sqrt(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) + \
           ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              (2*self.VY*(-(self.VX*self.Y) + self.X*self.VY) - \
                2*self.VZ*(self.VX*self.Z - self.X*self.VZ))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ))/\
            (2.*μ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              np.sqrt(1 - pow(-((self.VX*\
                      (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) + \
           (self.VX*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ))/\
            (μ*np.sqrt(1 - pow(-((self.VX*\
                       (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))))/\
         np.sqrt(1 - (pow(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ,2)*\
              (pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              pow(2/\
                 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ,2))/\
            (pow(μ,2)*(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              (pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2)))))/\
      np.sqrt(μ*pow(2/\
           np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
          (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,3)\
        )) - (3*μ*self.X*pow(2/\
         np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
        (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,2)*\
      (-(((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
             np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
               pow(self.VX*self.Z - self.X*self.VZ,2) + \
               pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
             (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
               (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/\
                μ))/\
           (μ*np.sqrt(1 - pow(-((self.VX*\
                      (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
               pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
               pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2)))) + \
        np.asin(((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
            np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
              pow(self.VX*self.Z - self.X*self.VZ,2) + \
              pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
            (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
              (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/\
               μ))/\
          (μ*np.sqrt(1 - pow(-((self.VX*\
                     (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) - \
              pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) - \
              pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2))*\
            np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                   μ) + self.X*\
                 (-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) + \
              pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) + \
              pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2))))))/\
    (pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5)*\
      pow(μ*pow(2/\
           np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
          (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,3)\
        ,1.5));

  aei_partials[6][2] = \
    -((((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             )*(-2*(-((self.VX*self.VY)/μ) + \
                (self.X*self.Y)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5))*(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                   μ) + self.X*\
                 (-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) - \
             2*(-(pow(self.VY,2)/μ) + \
                pow(self.Y,2)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - 1/\
                 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) + \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) - \
             2*((self.Y*self.Z)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - (self.VY*self.VZ)/μ)*\
              (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ))))/\
         (2.*μ*pow(1 - \
             pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2),1.5)) + \
        (2*self.Y*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
         (μ*pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5)*\
           np.sqrt(1 - pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                  μ) + self.X*\
                (-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                  (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) - \
        ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           (-2*self.VX*(-(self.VX*self.Y) + self.X*self.VY) + \
             2*self.VZ*(-(self.VY*self.Z) + self.Y*self.VZ))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ))/\
         (2.*μ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           np.sqrt(1 - pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                  μ) + self.X*\
                (-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                  (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) - \
        (self.VY*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ))/\
         (μ*np.sqrt(1 - pow(-((self.VX*\
                    (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) + \
        (-((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
               np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                 pow(self.VX*self.Z - self.X*self.VZ,2) + \
                 pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
               (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                 (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ)*\
               (2*(-((self.VX*self.VY)/μ) + \
                    (self.X*self.Y)/\
                     pow(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2),1.5))*\
                  (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) + \
                 2*(-(pow(self.VY,2)/μ) + \
                    pow(self.Y,2)/\
                     pow(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2),1.5) - \
                    1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                      + (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)*\
                  (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) + \
                 2*((self.Y*self.Z)/\
                     pow(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2),1.5) - (self.VY*self.VZ)/μ)*\
                  (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ))))/\
            (2.*μ*np.sqrt(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              pow(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2),1.5)) - \
           ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)*\
              (-2*(-((self.VX*self.VY)/μ) + \
                   (self.X*self.Y)/\
                    pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                     1.5))*\
                 (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) - \
                2*(-(pow(self.VY,2)/μ) + \
                   pow(self.Y,2)/\
                    pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                     1.5) - \
                   1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                    + (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)*\
                 (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) - \
                2*((self.Y*self.Z)/\
                    pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                     1.5) - (self.VY*self.VZ)/μ)*\
                 (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ))))/\
            (2.*μ*pow(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2),1.5)*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) - \
           (2*self.Y*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
            (μ*pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
               1.5)*np.sqrt(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) + \
           ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              (-2*self.VX*(-(self.VX*self.Y) + self.X*self.VY) + \
                2*self.VZ*(-(self.VY*self.Z) + self.Y*self.VZ))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ))/\
            (2.*μ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              np.sqrt(1 - pow(-((self.VX*\
                      (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) + \
           (self.VY*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ))/\
            (μ*np.sqrt(1 - pow(-((self.VX*\
                       (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))))/\
         np.sqrt(1 - (pow(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ,2)*\
              (pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              pow(2/\
                 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ,2))/\
            (pow(μ,2)*(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              (pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2)))))/\
      np.sqrt(μ*pow(2/\
           np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
          (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,3)\
        )) - (3*μ*self.Y*pow(2/\
         np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
        (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,2)*\
      (-(((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
             np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
               pow(self.VX*self.Z - self.X*self.VZ,2) + \
               pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
             (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
               (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/\
                μ))/\
           (μ*np.sqrt(1 - pow(-((self.VX*\
                      (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
               pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
               pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2)))) + \
        np.asin(((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
            np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
              pow(self.VX*self.Z - self.X*self.VZ,2) + \
              pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
            (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
              (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/\
               μ))/\
          (μ*np.sqrt(1 - pow(-((self.VX*\
                     (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) - \
              pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) - \
              pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2))*\
            np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                   μ) + self.X*\
                 (-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) + \
              pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) + \
              pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2))))))/\
    (pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5)*\
      pow(μ*pow(2/\
           np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
          (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,3)\
        ,1.5));

  aei_partials[6][3] =                                             \
    -((((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             )*(-2*((self.X*self.Z)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - (self.VX*self.VZ)/μ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) - \
             2*((self.Y*self.Z)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - (self.VY*self.VZ)/μ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) - \
             2*(pow(self.Z,2)/\
                 pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                  1.5) - 1/\
                 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                pow(self.VZ,2)/μ + \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)*\
              (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ))))/\
         (2.*μ*pow(1 - \
             pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2),1.5)) + \
        (2*self.Z*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
         (μ*pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5)*\
           np.sqrt(1 - pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                  μ) + self.X*\
                (-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                  (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) - \
        ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           (2*self.VX*(self.VX*self.Z - self.X*self.VZ) - \
             2*self.VY*(-(self.VY*self.Z) + self.Y*self.VZ))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ))/\
         (2.*μ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           np.sqrt(1 - pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                  μ) + self.X*\
                (-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                  (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) - \
        (self.VZ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ))/\
         (μ*np.sqrt(1 - pow(-((self.VX*\
                    (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) + \
        (-((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
               np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                 pow(self.VX*self.Z - self.X*self.VZ,2) + \
                 pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
               (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                 (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ)*\
               (2*((self.X*self.Z)/\
                     pow(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2),1.5) - (self.VX*self.VZ)/μ)*\
                  (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) + \
                 2*((self.Y*self.Z)/\
                     pow(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2),1.5) - (self.VY*self.VZ)/μ)*\
                  (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) + \
                 2*(pow(self.Z,2)/\
                     pow(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2),1.5) - \
                    1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                      - pow(self.VZ,2)/μ + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)*\
                  (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ))))/\
            (2.*μ*np.sqrt(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              pow(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2),1.5)) - \
           ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)*\
              (-2*((self.X*self.Z)/\
                    pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                     1.5) - (self.VX*self.VZ)/μ)*\
                 (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) - \
                2*((self.Y*self.Z)/\
                    pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                     1.5) - (self.VY*self.VZ)/μ)*\
                 (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) - \
                2*(pow(self.Z,2)/\
                    pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
                     1.5) - \
                   1/\
                    np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                    - pow(self.VZ,2)/μ + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)*\
                 (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ))))/\
            (2.*μ*pow(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2),1.5)*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) - \
           (2*self.Z*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
            (μ*pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),\
               1.5)*np.sqrt(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) + \
           ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              (2*self.VX*(self.VX*self.Z - self.X*self.VZ) - \
                2*self.VY*(-(self.VY*self.Z) + self.Y*self.VZ))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ))/\
            (2.*μ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              np.sqrt(1 - pow(-((self.VX*\
                      (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) + \
           (self.VZ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ))/\
            (μ*np.sqrt(1 - pow(-((self.VX*\
                       (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))))/\
         np.sqrt(1 - (pow(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ,2)*\
              (pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              pow(2/\
                 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ,2))/\
            (pow(μ,2)*(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              (pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2)))))/\
      np.sqrt(μ*pow(2/\
           np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
          (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,3)\
        )) - (3*μ*self.Z*pow(2/\
         np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
        (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,2)*\
      (-(((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
             np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
               pow(self.VX*self.Z - self.X*self.VZ,2) + \
               pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
             (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
               (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/\
                μ))/\
           (μ*np.sqrt(1 - pow(-((self.VX*\
                      (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
               pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
               pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2)))) + \
        np.asin(((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
            np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
              pow(self.VX*self.Z - self.X*self.VZ,2) + \
              pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
            (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
              (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/\
               μ))/\
          (μ*np.sqrt(1 - pow(-((self.VX*\
                     (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) - \
              pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) - \
              pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2))*\
            np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                   μ) + self.X*\
                 (-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) + \
              pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) + \
              pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2))))))/\
    (pow(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2),1.5)*\
      pow(μ*pow(2/\
           np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
          (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,3)\
        ,1.5));

  aei_partials[6][4] = \
    -((((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             )*(-2*((self.X*self.VX)/μ - \
                (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) - \
             2*((2*self.VX*self.Y)/μ - (self.X*self.VY)/μ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) - \
             2*((2*self.VX*self.Z)/μ - (self.X*self.VZ)/μ)*\
              (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ))))/\
         (2.*μ*pow(1 - \
             pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2),1.5)) + \
        (2*self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
         (pow(μ,2)*np.sqrt(1 - \
             pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) - \
        ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           (-2*self.Y*(-(self.VX*self.Y) + self.X*self.VY) + 2*self.Z*(self.VX*self.Z - self.X*self.VZ))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ))/\
         (2.*μ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           np.sqrt(1 - pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                  μ) + self.X*\
                (-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                  (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) - \
        (self.X*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ))/\
         (μ*np.sqrt(1 - pow(-((self.VX*\
                    (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) + \
        (-((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
               np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                 pow(self.VX*self.Z - self.X*self.VZ,2) + \
                 pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
               (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                 (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ)*\
               (2*((self.X*self.VX)/μ - \
                    (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
                  (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) + \
                 2*((2*self.VX*self.Y)/μ - (self.X*self.VY)/μ)*\
                  (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) + \
                 2*((2*self.VX*self.Z)/μ - (self.X*self.VZ)/μ)*\
                  (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ))))/\
            (2.*μ*np.sqrt(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              pow(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2),1.5)) - \
           ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)*\
              (-2*((self.X*self.VX)/μ - \
                   (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
                 (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) - \
                2*((2*self.VX*self.Y)/μ - (self.X*self.VY)/μ)*\
                 (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) - \
                2*((2*self.VX*self.Z)/μ - (self.X*self.VZ)/μ)*\
                 (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ))))/\
            (2.*μ*pow(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2),1.5)*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) - \
           (2*self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
            (pow(μ,2)*np.sqrt(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) + \
           ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              (-2*self.Y*(-(self.VX*self.Y) + self.X*self.VY) + \
                2*self.Z*(self.VX*self.Z - self.X*self.VZ))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ))/\
            (2.*μ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              np.sqrt(1 - pow(-((self.VX*\
                      (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) + \
           (self.X*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ))/\
            (μ*np.sqrt(1 - pow(-((self.VX*\
                       (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))))/\
         np.sqrt(1 - (pow(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ,2)*\
              (pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              pow(2/\
                 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ,2))/\
            (pow(μ,2)*(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              (pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2)))))/\
      np.sqrt(μ*pow(2/\
           np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
          (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,3)\
        )) - (3*self.VX*pow(2/\
         np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
        (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,2)*\
      (-(((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
             np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
               pow(self.VX*self.Z - self.X*self.VZ,2) + \
               pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
             (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
               (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/\
                μ))/\
           (μ*np.sqrt(1 - pow(-((self.VX*\
                      (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
               pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
               pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2)))) + \
        np.asin(((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
            np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
              pow(self.VX*self.Z - self.X*self.VZ,2) + \
              pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
            (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
              (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/\
               μ))/\
          (μ*np.sqrt(1 - pow(-((self.VX*\
                     (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) - \
              pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) - \
              pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2))*\
            np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                   μ) + self.X*\
                 (-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) + \
              pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) + \
              pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2))))))/\
    pow(μ*pow(2/\
         np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
        (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,3),\
     1.5);


  aei_partials[6][5] =  \
    -((((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             )*(-2*(-((self.VX*self.Y)/μ) + (2*self.X*self.VY)/μ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) - \
             2*((self.Y*self.VY)/μ - (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) - \
             2*((2*self.VY*self.Z)/μ - (self.Y*self.VZ)/μ)*\
              (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ))))/\
         (2.*μ*pow(1 - \
             pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2),1.5)) + \
        (2*self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
         (pow(μ,2)*np.sqrt(1 - \
             pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) - \
        ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           (2*self.X*(-(self.VX*self.Y) + self.X*self.VY) - \
             2*self.Z*(-(self.VY*self.Z) + self.Y*self.VZ))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ))/\
         (2.*μ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           np.sqrt(1 - pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                  μ) + self.X*\
                (-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                  (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) - \
        (self.Y*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ))/\
         (μ*np.sqrt(1 - pow(-((self.VX*\
                    (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) + \
        (-((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
               np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                 pow(self.VX*self.Z - self.X*self.VZ,2) + \
                 pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
               (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                 (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ)*\
               (2*(-((self.VX*self.Y)/μ) + (2*self.X*self.VY)/μ)*\
                  (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) + \
                 2*((self.Y*self.VY)/μ - \
                    (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
                  (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) + \
                 2*((2*self.VY*self.Z)/μ - (self.Y*self.VZ)/μ)*\
                  (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ))))/\
            (2.*μ*np.sqrt(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              pow(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2),1.5)) - \
           ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)*\
              (-2*(-((self.VX*self.Y)/μ) + (2*self.X*self.VY)/μ)*\
                 (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) - \
                2*((self.Y*self.VY)/μ - \
                   (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
                 (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) - \
                2*((2*self.VY*self.Z)/μ - (self.Y*self.VZ)/μ)*\
                 (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ))))/\
            (2.*μ*pow(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2),1.5)*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) - \
           (2*self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
            (pow(μ,2)*np.sqrt(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) + \
           ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              (2*self.X*(-(self.VX*self.Y) + self.X*self.VY) - \
                2*self.Z*(-(self.VY*self.Z) + self.Y*self.VZ))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ))/\
            (2.*μ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              np.sqrt(1 - pow(-((self.VX*\
                      (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) + \
           (self.Y*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ))/\
            (μ*np.sqrt(1 - pow(-((self.VX*\
                       (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))))/\
         np.sqrt(1 - (pow(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ,2)*\
              (pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              pow(2/\
                 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ,2))/\
            (pow(μ,2)*(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              (pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2)))))/\
      np.sqrt(μ*pow(2/\
           np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
          (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,3)\
        )) - (3*self.VY*pow(2/\
         np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
        (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,2)*\
      (-(((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
             np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
               pow(self.VX*self.Z - self.X*self.VZ,2) + \
               pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
             (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
               (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/\
                μ))/\
           (μ*np.sqrt(1 - pow(-((self.VX*\
                      (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
               pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
               pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2)))) + \
        np.asin(((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
            np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
              pow(self.VX*self.Z - self.X*self.VZ,2) + \
              pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
            (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
              (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/\
               μ))/\
          (μ*np.sqrt(1 - pow(-((self.VX*\
                     (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) - \
              pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) - \
              pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2))*\
            np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                   μ) + self.X*\
                 (-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) + \
              pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) + \
              pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2))))))/\
    pow(μ*pow(2/\
         np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
        (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,3),\
     1.5);

  aei_partials[6][6] =                                              \
    -((((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             )*(-2*(-((self.VX*self.Z)/μ) + (2*self.X*self.VZ)/μ)*\
              (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) - \
             2*(-((self.VY*self.Z)/μ) + (2*self.Y*self.VZ)/μ)*\
              (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ)) - \
             2*((self.Z*self.VZ)/μ - (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
              (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ))))/\
         (2.*μ*pow(1 - \
             pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2),1.5)) + \
        (2*self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
         (pow(μ,2)*np.sqrt(1 - \
             pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) - \
        ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
           (-2*self.X*(self.VX*self.Z - self.X*self.VZ) + 2*self.Y*(-(self.VY*self.Z) + self.Y*self.VZ))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ))/\
         (2.*μ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           np.sqrt(1 - pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                  μ) + self.X*\
                (-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                  (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) - \
        (self.Z*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
             pow(self.VX*self.Z - self.X*self.VZ,2) + \
             pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
           (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
             (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ\
             ))/\
         (μ*np.sqrt(1 - pow(-((self.VX*\
                    (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.X*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Y*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2) - \
             pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
               self.Z*(-(1/\
                     np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2))\
                     ) + (pow(self.VX,2) + pow(self.VY,2) + \
                     pow(self.VZ,2))/μ),2))) + \
        (-((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
               np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                 pow(self.VX*self.Z - self.X*self.VZ,2) + \
                 pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
               (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                 (pow(self.VX,2) + pow(self.VY,2) + \
                    pow(self.VZ,2))/μ)*\
               (2*(-((self.VX*self.Z)/μ) + (2*self.X*self.VZ)/μ)*\
                  (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) + \
                 2*(-((self.VY*self.Z)/μ) + (2*self.Y*self.VZ)/μ)*\
                  (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) + \
                 2*((self.Z*self.VZ)/μ - \
                    (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
                  (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                    self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                       (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ))))/\
            (2.*μ*np.sqrt(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              pow(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2),1.5)) - \
           ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ)*\
              (-2*(-((self.VX*self.Z)/μ) + (2*self.X*self.VZ)/μ)*\
                 (-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) - \
                2*(-((self.VY*self.Z)/μ) + (2*self.Y*self.VZ)/μ)*\
                 (-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ)) - \
                2*((self.Z*self.VZ)/μ - \
                   (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)/μ)*\
                 (-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                   self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                      (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ))))/\
            (2.*μ*pow(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2),1.5)*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) - \
           (2*self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2)))/\
            (pow(μ,2)*np.sqrt(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) + \
           ((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
              (-2*self.X*(self.VX*self.Z - self.X*self.VZ) + \
                2*self.Y*(-(self.VY*self.Z) + self.Y*self.VZ))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ))/\
            (2.*μ*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              np.sqrt(1 - pow(-((self.VX*\
                      (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))) + \
           (self.Z*np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ))/\
            (μ*np.sqrt(1 - pow(-((self.VX*\
                       (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))))/\
         np.sqrt(1 - (pow(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ,2)*\
              (pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
                pow(self.VX*self.Z - self.X*self.VZ,2) + \
                pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
              pow(2/\
                 np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
                (pow(self.VX,2) + pow(self.VY,2) + \
                   pow(self.VZ,2))/μ,2))/\
            (pow(μ,2)*(1 - \
                pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2))*\
              (pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                  self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) + \
                pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                     μ) + \
                  self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                     (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2)))))/\
      np.sqrt(μ*pow(2/\
           np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
          (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,3)\
        )) - (3*self.VZ*pow(2/\
         np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
        (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,2)*\
      (-(((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
             np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
               pow(self.VX*self.Z - self.X*self.VZ,2) + \
               pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
             (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
               (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/\
                μ))/\
           (μ*np.sqrt(1 - pow(-((self.VX*\
                      (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.X*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
               pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.Y*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2) - \
               pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                 self.Z*(-(1/\
                       np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                       pow(self.Z,2))) + \
                    (pow(self.VX,2) + pow(self.VY,2) + \
                       pow(self.VZ,2))/μ),2)))) + \
        np.asin(((self.X*self.VX + self.Y*self.VY + self.Z*self.VZ)*\
            np.sqrt(pow(-(self.VX*self.Y) + self.X*self.VY,2) + \
              pow(self.VX*self.Z - self.X*self.VZ,2) + \
              pow(-(self.VY*self.Z) + self.Y*self.VZ,2))*\
            (2/np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
              (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/\
               μ))/\
          (μ*np.sqrt(1 - pow(-((self.VX*\
                     (self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.X*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) - \
              pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) - \
              pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2))*\
            np.sqrt(pow(-((self.VX*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/\
                   μ) + self.X*\
                 (-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) + \
              pow(-((self.VY*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Y*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2) + \
              pow(-((self.VZ*(self.X*self.VX + self.Y*self.VY + self.Z*self.VZ))/μ) + \
                self.Z*(-(1/\
                      np.sqrt(pow(self.X,2) + pow(self.Y,2) + \
                      pow(self.Z,2))) + \
                   (pow(self.VX,2) + pow(self.VY,2) + \
                      pow(self.VZ,2))/μ),2))))))/\
    pow(μ*pow(2/\
         np.sqrt(pow(self.X,2) + pow(self.Y,2) + pow(self.Z,2)) - \
        (pow(self.VX,2) + pow(self.VY,2) + pow(self.VZ,2))/μ,3),\
     1.5);

  return aei_partials
