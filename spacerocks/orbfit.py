import numpy as np

class OrbitFit:

    def __init__(self, detections):
        pass
        #if mpc_format == True:
        #    self.file = file

        #elif spacerocks_format = True:


    #def equatorial_to_tangent(self):
    #    '''
    #    Project the equatorial position to the tangent plane
    #    '''
    #    β = np.arcsin(np.cos(ε) * np.sin(self.dec) -  \
    #        np.sin(ε) * np.cos(self.dec) * np.sin(self.ra))
    #    λ = np.arctan((np.cos(ε) * np.cos(self.dec) * np.sin(self.ra) + \
    #        np.sin(ε) * np.sin(self.dec)) / (np.cos(self.dec) * np.cos(self.ra)))


    #    #dβ_ddec = (np.cos(self.dec)*np.cos(ε) + np.sin(self.ra)*np.sin(self.dec)*np.sin(ε))**2 / \
    #    #          (1 - (np.cos(ε)*np.sin(self.dec) - np.cos(self.dec)*np.sin(self.ra)*np.sin(ε))**2)

    #    #dβ_dra = (np.cos(self.ra) * np.cos(self.dec) * np.sin(ε))**2 / \
    #    #         (1 - (np.cos(ε)*np.sin(self.dec) - np.cos(self.dec)*np.sin(self.ra)*np.sin(ε))**2)

    #    #dλ_ddec = (np.cos(self.ra)**2 * np.sec(self.dec)**4 * np.sin(ε)**2) \
    #    #          / (np.cos(self.ra)**2 + (np.cos(ε)*np.sin(self.ra) + np.sin(ε)*np.tan(self.dec))**2)**2

    #    #dλ_dra = (np.cos(ε) + np.sin(self.ra) * np.sin(ε) * np.tan(self.dec))**2 \
    #    #         / (np.cos(self.ra)**2 + np.cos(ε)**2 * np.sin(self.ra)**2) \
    #    #         + np.tan(self.dec) * (np.sin(self.ra) * np.sin(2*ε) \
    #    #         + np.sin(ε)**2 * np.tan(self.dec)))**2

    #    return θx, θy#, σθx, σθy

    @property
    def beta(self):
        return arcsin(cos(obliquity) * sin(dec) - sin(obliquity) * cos(dec) * sin(ra))

    @property
    def theta_x(self):
        #b =
        #l =
        return cos(b) * sin(l - l[0])) / (sin(b[0]) * sin(b) - cos(b[0]) * cos(b) * cos(l - l[0]))

    @property
    def theta_y(self):
        #b =
        #l =
        return cos(b[0]) * cos(b) - sin(b[0]) * np.sin(b) * cos(l - l[0]))/ (sin(b[0]) * sin(b) + cos(b[0]) * cos(b) * np.cos(l - l[0]))
