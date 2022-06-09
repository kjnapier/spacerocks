from .units import Units

import numpy as np
import dateutil
import datetime

from astropy.time import Time
from astropy.coordinates import Angle

def great_circle_distance(ra1: Angle, dec1: Angle, ra2: Angle, dec2: Angle):
    return np.arccos(np.sin(dec1.rad) * np.sin(dec2.rad) + 
                     np.cos(dec1.rad) * np.cos(dec2.rad) * np.cos(ra1.rad - ra2.rad))

def time_handler(d, units=Units()):
    d = np.atleast_1d(d)
    if units.timeformat is not None:
        epoch = Time(d, format=units.timeformat, scale=units.timescale)
    else:
        epoch = infer_time_format(d, units)
    
    return epoch

def infer_time_format(d, units=Units()):
    if isinstance(d[0], Time):
        epoch = d
    elif isinstance(d[0], datetime.datetime):
        epoch = Time(d, format='datetime', scale=units.timescale)
    elif isinstance(d[0], str):
        dates = [dateutil.parser.parse(x, fuzzy_with_tokens=True)[0] for x in d]
        epoch = Time(dates, format='datetime', scale=units.timescale)  
    elif np.isscalar(d[0]):
        if np.all(d > 100000):
            epoch = Time(d, format='jd', scale=units.timescale)
        else:
            epoch = Time(d, format='mjd', scale=units.timescale)
            
    return epoch

