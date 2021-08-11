from astropy.coordinates import Angle
from astropy import units as u

class Prediction:

    def __init__(self, **kwargs):

        self.epoch = kwargs.get('epoch')
        self.ra = Angle(kwargs.get('ra'), u.rad)
        self.dec = Angle(kwargs.get('dec'), u.rad)

        self.err_a = Angle(kwargs.get('err_a'), u.arcsec)
        self.err_b = Angle(kwargs.get('err_b'), u.arcsec)
        self.err_pa = Angle(kwargs.get('err_pa'), u.deg)

        self.elong = Angle(kwargs.get('elong'), u.deg)
        self.opp = Angle(kwargs.get('opp'), u.deg)
