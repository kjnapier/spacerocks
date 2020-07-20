class OrbitFit:

    def __init__(self, ra, dec, epoch, obscode):

        self.ra = ra
        self.dec = dec

        pass

    def equatorial_to_tangent(self):
        '''
        Project the equatorial position to the tangent plane
        '''
        β = np.arcsin(np.cos(ε) * np.sin(self.dec) -  \
            np.sin(ε) * np.cos(self.dec) * np.sin(self.ra))
        λ = np.arctan((np.cos(ε) * np.cos(self.dec) * np.sin(self.ra) + \
            np.sin(ε) * np.sin(self.dec)) / (np.cos(self.dec) * np.cos(self.ra)))

        θx = (np.cos(β) * np.sin(λ - λ[0])) / (np.sin(b[0]) * np.sin(β) - \
                np.cos(β[0]) * np.cos(β) * np.cos(λ - λ[0]))
        θy = (np.cos(β[0]) * np.cos(β) - np.sin(β[0]) * np.sin(β) * np.cos(λ - λ[0]))/ \
             (np.sin(β[0]) * np.sin(β) + np.cos(β[0]) * np.cos(β) * np.cos(λ - λ[0]))

        #dβ_ddec = (np.cos(self.dec)*np.cos(ε) + np.sin(self.ra)*np.sin(self.dec)*np.sin(ε))**2 / \
        #          (1 - (np.cos(ε)*np.sin(self.dec) - np.cos(self.dec)*np.sin(self.ra)*np.sin(ε))**2)

        #dβ_dra = (np.cos(self.ra) * np.cos(self.dec) * np.sin(ε))**2 / \
        #         (1 - (np.cos(ε)*np.sin(self.dec) - np.cos(self.dec)*np.sin(self.ra)*np.sin(ε))**2)

        #dλ_ddec = (np.cos(self.ra)**2 * np.sec(self.dec)**4 * np.sin(ε)**2) \
        #          / (np.cos(self.ra)**2 + (np.cos(ε)*np.sin(self.ra) + np.sin(ε)*np.tan(self.dec))**2)**2

        #dλ_dra = (np.cos(ε) + np.sin(self.ra) * np.sin(ε) * np.tan(self.dec))**2 \
        #         / (np.cos(self.ra)**2 + np.cos(ε)**2 * np.sin(self.ra)**2) \
        #         + np.tan(self.dec) * (np.sin(self.ra) * np.sin(2*ε) \
        #         + np.sin(ε)**2 * np.tan(self.dec)))**2

        return θx, θy#, σθx, σθy
