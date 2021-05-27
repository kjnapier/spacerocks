#import numpy as np
#
#class OrbitFit:
#
#    def __init__(self, detections):
#        pass
#        #if mpc_format == True:
#        #    self.file = file
#
#        #elif spacerocks_format = True:
#
#    '''
#    We'll have 4-5 well-measured quantities:
#    alpha, beta, gamma, alpha_dot, and beta_dot.
#    gamma_dot is unknown, but there is an energy constraint.
#
#    We have well-defined uncertainties for ra, dec, ra rate, and dec rate.
#
#    How do we get uncertainties on distance?
#    I think those come from ra rate and dec rate. Maybe we can write down the
#    distance equation, and then propagate the uncertainties.
#
#
#    '''
#
#    @property
#    def beta(self):
#        return arcsin(cos(obliquity) * sin(dec) - sin(obliquity) * cos(dec) * sin(ra))
#
#    @property
#    def lambda(self):
#        return arctan((cos(obliquity) * cos(dec) * sin(ra) + sin(obliquity) * sin(dec)) / (cos(dec) * cos(ra)))
#
#    @property
#    def theta_x(self):
#        #b =
#        #l =
#        return cos(b) * sin(l - l[0])) / (sin(b[0]) * sin(b) - cos(b[0]) * cos(b) * cos(l - l[0]))
#
#    @property
#    def theta_y(self):
#        #b =
#        #l =
#        return cos(b[0]) * cos(b) - sin(b[0]) * np.sin(b) * cos(l - l[0]))/ (sin(b[0]) * sin(b) + cos(b[0]) * cos(b) * np.cos(l - l[0]))
