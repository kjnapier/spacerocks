import os
import warnings

import healpy as hp
import numpy as np
import pandas as pd
from skyfield.api import Topos, Loader
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord


from .constants import epsilon, mu_bary, c
from .transformations import Transformations
from .linalg3d import norm
from .convenience import Convenience


try:
    import cartopy.crs as ccrs
except:
    pass

# Read in the observatory codes file and rehash as a dataframe.
observatories = pd.read_csv(os.path.join(os.path.dirname(__file__),
                            'data',
                            'observatories.csv'))

# Load in planets for ephemeride calculation.
load = Loader('./Skyfield-Data', expire=False, verbose=False)
ts = load.timescale()
planets = load('de423.bsp')

class Observe(Transformations, Convenience):

    def __init__(self, rocks, obscode=None, NSIDE=None):

        if obscode is not None:
            Observe.obscode = str(obscode).zfill(3)
            obs = observatories[observatories.obscode == Observe.obscode]
            Observe.obslat = obs.lat.values
            Observe.obslon = obs.lon.values
            Observe.obselev = obs.elevation.values
        else:
            Observe.obscode = 'Geocenter'

        if NSIDE is not None:
            if np.isscalar(NSIDE):
                NSIDE = np.array([NSIDE])
            Observe.NSIDE = NSIDE
        else:
            Observe.NSIDE = None

        self.xyz_to_equa(rocks)

        if rocks.frame == 'barycentric':
            self.xyz_to_equa(rocks)

        elif rocks.frame == 'heliocentric':
            rocks.to_bary()
            self.xyz_to_equa(rocks)
            rocks.to_helio()

        try:
            self.mag = self.estimate_mag(rocks)
        except:
            pass

        if Observe.NSIDE is not None:
            for value in Observe.NSIDE:
                setattr(self, 'HPIX_{}'.format(value), self.radec_to_hpix(value))

        self.name = rocks.name
        self.epoch = rocks.epoch


    def xyz_to_equa(self, rocks):
        '''
        Transform from barycentric Cartesian coordinates to equatorial
        coordinates. If you need very precise values, the method will use
        a topocentric correction to the Earth's position.
        See https://en.wikipedia.org/wiki/Horizontal_coordinate_system.
        This results in a significant slowdown to the code.
        '''
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            t = ts.tt(jd=rocks.epoch.tt.jd)
            earth = planets['earth']

            # Only used for the topocentric calculation.
            if Observe.obscode != 'Geocenter':
                earth += Topos(latitude_degrees=Observe.obslat,
                                longitude_degrees=Observe.obslon,
                                elevation_m=Observe.obselev) # topocentric calculation

            x_earth, y_earth, z_earth = earth.at(t).position.au * u.au # earth ICRS position
            earth_dis = norm([x_earth, y_earth, z_earth])
            x0, y0, z0 = rocks.x, rocks.y, rocks.z
            for idx in range(10):
                # transfer ecliptic to ICRS and shift to Geocentric (topocentric)
                x = x0 - x_earth
                y = y0 * np.cos(epsilon) - z0 * np.sin(epsilon) - y_earth
                z = y0 * np.sin(epsilon) + z0 * np.cos(epsilon) - z_earth
                delta = norm([x, y, z])
                ltt = delta / c
                M = rocks.M - ltt * (mu_bary / rocks.a**3)**0.5
                if idx < 9:
                    x0, y0, z0 = self.kep_to_xyz_pos(rocks.a, rocks.e, rocks.inc,
                                                    rocks.arg, rocks.node, M)

            # Cartesian to spherical coordinate
            self.delta = norm([x, y, z])
            self.ltt = self.delta / c
            self.dec = Angle(np.arcsin(z / norm([x, y, z])), u.rad)
            self.ra = Angle(np.arctan2(y, x), u.rad).wrap_at(2 * np.pi * u.rad)
            self.phase_angle = Angle(np.arccos(-(earth_dis**2 - rocks.r**2 - self.delta**2)/(2 * rocks.r* self.delta)), u.rad)
            self.elong = Angle(np.arccos(-(rocks.r**2 - self.delta**2 - earth_dis**2)/(2 * self.delta * earth_dis)), u.rad)
            self.skycoord = SkyCoord(self.ra, self.dec, frame='icrs')

        return self

    def radec_to_hpix(self, NSIDE):
        '''
        Convert (ra, dec) into healpix.
        '''
        return hp.pixelfunc.ang2pix(NSIDE, np.pi/2 - self.dec.radian, self.ra.radian, nest=True)

    def estimate_mag(self, rocks):
        '''
        Estimate the apparent magnitude of a TNO
        https://iopscience.iop.org/article/10.3847/1538-3881/ab18a9/pdf for light curves
        '''
        q = (rocks.r**2 + self.delta**2 - 1 * u.au**2)/(2 * rocks.r * self.delta)

        ## pyephem
        beta = np.zeros(len(q))
        beta[np.where(q <= 1)[0]] = np.pi * u.rad
        beta[np.where(q >= 1)[0]] = 0 * u.rad

        Psi_1 = np.exp(-3.33 * np.tan(beta/2)**0.63)
        Psi_2 = np.exp(-1.87 * np.tan(beta/2)**1.22)
        mag = rocks.H + 5 * np.log10(rocks.r * self.delta / u.au**2)

        not_zero = np.where((Psi_1 != 0) | (Psi_2 != 0))[0]
        mag[not_zero] -= 2.5 * np.log10((1 - rocks.G[not_zero]) * Psi_1[not_zero] + rocks.G[not_zero] * Psi_2[not_zero])

        #try:
        #    mag += rocks.light_amp * np.sin( + rocks.phi0)
        #except:
        #    pass

        return mag

    def plot_radec(self, color='black', alpha=0.5, zoom=False, galactic_plane=False, ecliptic_plane=True):
        '''
        Plot the right ascension and declination of each object on a Mollweide
        projection. (See https://en.wikipedia.org/wiki/Mollweide_projection)
        '''
        if zoom == True:

            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

            xdata = self.ra.degree
            xdata[xdata > 180] -= 360

            xmin=np.min(xdata)
            xmax=np.max(xdata)
            ymin=np.min(self.dec.degree)
            ymax=np.max(self.dec.degree)

            ax.scatter(-xdata, self.dec.degree, color='black', alpha=0.5)

            xticks = np.linspace(-xmax, -xmin, 8)
            yticks = np.linspace(ymin, ymax, 8)

            ax.set_xticks(xticks)
            ax.set_yticks(yticks)
            xticklabels = [r'${:.2f}\degree$'.format(-value) for value in xticks]
            yticklabels = [r'${:.2f}\degree$'.format(value) for value in yticks]
            ax.set_xticklabels(xticklabels)
            ax.set_yticklabels(yticklabels)

            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                              linewidth=1, color='gray', alpha=0.5, linestyle='-')
            gl.xlocator = mticker.FixedLocator(xticks)
            gl.ylocator = mticker.FixedLocator(yticks)

            gl.bottom_labels = False
            gl.top_labels = False
            gl.left_labels = False
            gl.right_labels = False

            xrange = xmax - xmin
            yrange = ymax - ymin

            try:
                ax.set_extent([-xmax - xrange * 0.05, -xmin + xrange * 0.05,
                               ymin - yrange * 0.05, ymax + yrange * 0.05], crs=ccrs.PlateCarree())
            except:
                ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

            if ecliptic_plane == True:
                def radec2project(ra, dec):
                    ra[ra>180] -= 360
                    return (ra, dec)

                plane_lon = np.linspace(-np.pi, np.pi, 1000)
                plane_lat = np.zeros(len(plane_lon))
                ecl_plane = np.zeros([len(plane_lat), 2])

                for i in range(len(plane_lat)):
                    ecl_plane[i] = ephem.Equatorial(ephem.Ecliptic(plane_lon[i], plane_lat[i])).get()

                x, y = radec2project(np.degrees(ecl_plane.T[0]), np.degrees(ecl_plane.T[1]))
                ax.plot(x[1:999], -y[1:999], 'r-', zorder=2)

            if galactic_plane == True:

                galactic = SkyCoord(l=np.linspace(0, 360, 1000)*u.degree, b=np.zeros(1000)*u.degree, frame='galactic')
                galxdata = galactic.icrs.ra.degree
                galxdata[galxdata > 180] -= 360
                galydata = galactic.icrs.dec.degree
                order = galxdata.argsort()
                ax.plot(-galxdata[order][1:999], -galydata[order][1:999], 'g-', zorder=3)

        else:

            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111, projection='mollweide')
            ax.grid(True)

            xdata = self.ra
            xdata[xdata.value > np.pi] -= 2*np.pi * u.rad
            ax.set_xticklabels([r'$150\degree$', r'$120\degree$', r'$90\degree$' ,
                                r'$60\degree$' , r'$30\degree$' , r'$0\degree$'  ,
                                r'$330\degree$', r'$300\degree$', r'$270\degree$',
                                r'$240\degree$', r'$210\degree$'])

            ax.scatter(-xdata, self.dec, color=color, alpha=alpha)

            if ecliptic_plane == True:
                def radec2project(ra, dec):
                    ra[ra>180] -= 360
                    return (ra, dec)

                plane_lon = np.linspace(-np.pi, np.pi, 1000)
                plane_lat = np.zeros(len(plane_lon))
                ecl_plane = np.zeros([len(plane_lat), 2])

                for i in range(len(plane_lat)):
                    ecl_plane[i] = ephem.Equatorial(ephem.Ecliptic(plane_lon[i], plane_lat[i])).get()

                x, y = radec2project(np.degrees(ecl_plane.T[0]), np.degrees(ecl_plane.T[1]))
                ax.plot(np.radians(x[1:999]), -np.radians(y[1:999]), 'r-', zorder=2)

            if galactic_plane == True:

                galactic = SkyCoord(l=np.linspace(0, 360, 1000)*u.degree, b=np.zeros(1000)*u.degree, frame='galactic')
                galxdata = galactic.icrs.ra.rad
                galxdata[galxdata > np.pi] -= 2*np.pi
                galydata = galactic.icrs.dec.rad
                order = galxdata.argsort()
                ax.plot(-galxdata[order][1:999], galydata[order][1:999], 'g-', zorder=3)

        return fig, ax
