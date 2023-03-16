from .image import Image
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import Angle, SkyCoord
from astropy import units as u
from astropy.time import Time
import numpy as np

class ExposureSequence:

    def __init__(self, diffim_files, weight_files, hardware='CPU'):
       
        self.hardware = hardware
        self.diffims = [Image(x, hardware=self.hardware) for x in diffim_files]
        self.weights = [Image(x, hardware=self.hardware) for x in weight_files]   
        self.epoch = Time(np.array([im.epoch.jd for im in self.diffims]), format='jd')
        self.fwhm = [im.fwhm for im in self.diffims]