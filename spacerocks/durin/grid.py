import numpy as np

from astropy.coordinates import Angle
from astropy.time import Time
from astropy import units as u

from ..constants import mu_bary
from ..spice import SpiceBody
from ..units import Units
from ..spacerock import SpaceRock

import alphashape

def compute_hull(ra: Angle, dec: Angle, epoch: Time, rmin: float = 30, rmax: float = 1000, imin: float = 0, imax: float = 180) -> alphashape.alphashape:

    mu = mu_bary.value
    earth = SpiceBody(spiceid='Earth')
    e = earth.at(epoch).position
    x_earth, y_earth, z_earth = e.x.au[0], e.y.au[0], e.z.au[0]
    
    N = 10_000
    phi = np.random.uniform(-np.pi, np.pi, N)
    theta = np.arccos(2 * np.random.uniform(0, 1, N) - 1.0)
    
    ra_rates = []
    dec_rates = []
    
    units = Units()
    units.timescale = 'utc'
    units.timeformat = 'jd'
    
    rmin = rmin
    rmax = rmax
    
    imin = imin
    imax = imax
    
    for r in np.logspace(np.log10(rmin), np.log10(rmax), 50):
    
        x = (r * np.cos(dec.rad) * np.cos(ra.rad) + x_earth)
        y = (r * np.cos(dec.rad) * np.sin(ra.rad) + y_earth)
        z = (r * np.sin(dec.rad) + z_earth)
        
        rbary = np.sqrt(x*x + y*y + z*z)
        
        vmax = np.sqrt(mu * (2 / rbary - 1/2_000_000))
        v = vmax * np.random.rand(N)**(1/3)
        vx = v * np.cos(theta) * np.sin(phi)
        vy = v * np.sin(theta) * np.sin(phi)
        vz = v * np.cos(phi)
        
        vx = vx.flatten()
        vy = vy.flatten()
        vz = vz.flatten()
        
        N = len(vx)
        
        rr = SpaceRock(x=np.repeat(x, N), 
                       y=np.repeat(y, N), 
                       z=np.repeat(z, N), 
                       vx=vx, 
                       vy=vy, 
                       vz=vz, 
                       epoch=np.repeat(epoch.utc.jd, N), 
                       name=np.array([idx for idx in range(N)], dtype='object'),
                       frame='J2000',
                       origin='ssb', 
                       units=units)
        
        rr.change_frame('ECLIPJ2000')
        rr = rr[(rr.inc.deg > imin) * (rr.inc.deg < imax)]
                
        oo = rr.observe(spiceid='Earth')
        ra_rates.extend(oo.ra_rate.to(u.arcsec/u.hour).value)
        dec_rates.extend(oo.dec_rate.to(u.arcsec/u.hour).value)
    
    ra_rates = np.ndarray.flatten(np.array(ra_rates))
    dec_rates = np.ndarray.flatten(np.array(dec_rates))
    points = np.array([ra_rates, dec_rates]).T
    hull = alphashape.alphashape(points, 0.1)

    return hull

def compute_stack_rates(hull: alphashape.alphashape, maximum_trailing: float, delta_t: float):
    
    step = (maximum_trailing * u.arcsec  / delta_t).to(u.arcsec / u.hour).value

    xmin, ymin, xmax, ymax = hull.bounds
    xmin -= 1
    ymin -= 1
    xmax += 1
    ymax += 1

    X = np.arange(xmin, xmax, step)
    Y = np.arange(ymin, ymax, step)
    x, y = np.meshgrid(X, Y)
    x[::2] += np.diff(X)[0]/2
    x = x.ravel()
    y = y.ravel()

    in_hull = []
    for xx, yy in zip(x, y):
        if hull.contains(alphashape.alphashape(np.array([[xx, yy]]))):
            in_hull.append(True)
        else:
            in_hull.append(False)

    ra_rates = x[in_hull]
    dec_rates = y[in_hull]

    return ra_rates, dec_rates