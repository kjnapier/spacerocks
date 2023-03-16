from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import Angle, SkyCoord
from astropy import units as u
from astropy.time import Time

class Image:
    '''
    Helper class for working with the DEEP fits images.

    '''
    def __init__(self, file, hardware='CPU'):
        hdu = fits.open(file, memmap=False)
        if hardware == 'GPU':
            import cupy as cp
            if len(hdu) == 2:
                self.data = cp.array(hdu[1].data.copy())
                self.header = hdu[1].header
            else:
                self.data = cp.array(hdu[0].data.copy())
                self.header = hdu[0].header
        elif hardware == 'CPU':
            if len(hdu) == 2:
                self.data = hdu[1].data.copy()
                self.header = hdu[1].header
            else:
                self.data = hdu[0].data.copy()
                self.header = hdu[0].header

        self.__calc_ra_dec()
        
    @property
    def wcs(self):
        if not hasattr(self, '_wcs'):
            self.wcs = WCS(self.header)
        return self._wcs

    @wcs.setter
    def wcs(self, value):
        self._wcs = value

    @wcs.deleter
    def wcs(self):
        del self._wcs

    @property
    def fwhm(self):
        return 4 #self.header['FWHM']

    @property
    def epoch(self):
        return Time(self.header['MJD-OBS'], format='mjd', scale='utc')

    @property
    def expnum(self):
        return self.header['EXPNUM']

    def __calc_ra_dec(self):
        Y_size, X_size = self.data.shape
        Ycenter = round(Y_size/2)
        Xcenter = round(X_size/2)
        ra, dec = self.wcs.wcs_pix2world(Xcenter, Ycenter, 0)
        self.ra = Angle(ra, u.deg)
        self.dec = Angle(dec, u.deg)