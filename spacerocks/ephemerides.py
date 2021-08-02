import random
import copy
import os

from astropy import units as u
from astropy.coordinates import Angle, Distance
from astropy.time import Time

from numpy import sin, cos, arctan2, sqrt, array, pi, zeros, arcsin, tan
import numpy as np
import pandas as pd

from .constants import c, epsilon

from skyfield.api import Topos, Loader
# Load in planets for ephemeride calculation.
load = Loader('./Skyfield-Data', expire=False, verbose=False)
ts = load.timescale()
planets = load('de440.bsp')
earth = planets['earth']
sun = planets['sun']

from .convenience import Convenience

class Ephemerides(Convenience):

    '''
    Take in the ecliptic state vector (relative to the observer) and
    return any quantity
    '''

    def __init__(self, **kwargs):
        self.x = kwargs.get('x')
        self.y = kwargs.get('y')
        self.z = kwargs.get('z')
        self.vx = kwargs.get('vx')
        self.vy = kwargs.get('vy')
        self.vz = kwargs.get('vz')
        self.epoch = kwargs.get('epoch')
        self.name = kwargs.get('name')

        if kwargs.get('H0') is not None:
            self.H0 = kwargs.get('H0')
            self.G = kwargs.get('G')

        if kwargs.get('r_helio') is not None:
            self.r_helio = kwargs.get('r_helio')

        if kwargs.get('delta_H') is not None:
            self.delta_H = kwargs.get('delta_H')
            self.rotation_period = kwargs.get('rotation_period')
            self.phi0 = kwargs.get('phi0')
            self.t0 = kwargs.get('t0')

        self.delta = Distance(np.sqrt(self.x**2 + self.y**2 + self.z**2), u.au)



    @property
    def ra(self):
        if not hasattr(self, '_ra'):
            self.ra = Angle(np.arctan2(self.y, self.x), u.rad).wrap_at(2 * np.pi * u.rad)
        return self._ra

    @ra.setter
    def ra(self, value):
        self._ra = value

    @ra.deleter
    def ra(self):
        del self._ra

    @property
    def dec(self):
        if not hasattr(self, '_dec'):
            self.dec = Angle(np.arcsin(self.z / sqrt(self.x**2 + self.y**2 + self.z**2)), u.rad)
        return self._dec

    @dec.setter
    def dec(self, value):
        self._dec = value

    @dec.deleter
    def dec(self):
        del self._dec

    @property
    def ra_rate(self):
        if not hasattr(self, '_ra_rate'):
            self.ra_rate = -cos(self.dec) * (self.y * self.vx - self.x * self.vy) / (self.x**2 + self.y**2) * u.rad
        return self._ra_rate

    @ra_rate.setter
    def ra_rate(self, value):
        self._ra_rate = value

    @ra_rate.deleter
    def ra_rate(self):
        del self._ra_rate

    @property
    def dec_rate(self):
        if not hasattr(self, '_dec_rate'):
            self.dec_rate = (-self.z * (self.x * self.vx + self.y * self.vy) + ((self.x**2 + self.y**2) * self.vz)) \
                            / (sqrt(self.x**2 + self.y**2) * (self.x**2 + self.y**2 + self.z**2)) * u.rad
        return self._dec_rate

    @dec_rate.setter
    def dec_rate(self, value):
        self._dec_rate = value

    @dec_rate.deleter
    def dec_rate(self):
        del self._dec_rate

    @property
    def mag(self):
        if not hasattr(self, '_mag'):
            self.mag = self.estimate_mag()
        return self._mag

    @mag.setter
    def mag(self, value):
        self._mag = value

    @mag.deleter
    def mag(self):
        del self._mag


    @property
    def H(self):
        if not hasattr(self, '_H'):
            if hasattr(self, 'delta_H'):
                ltt = self.delta / c
                dH = self.delta_H * np.sin((self.epoch.jd - self.t0.jd - ltt.value) * 2 * np.pi / self.rotation_period + self.phi0.rad)

                self.H = self.H0 + dH

            else:
                self.H = self.H0
        return self._H

    @H.setter
    def H(self, value):
        self._H = value

    @H.deleter
    def H(self):
        del self._H


    @property
    def b(self):
        return Angle(arcsin(cos(epsilon) * sin(self.dec) - sin(epsilon) * cos(self.dec) * sin(self.ra)), u.rad)

    @property
    def l(self):
        return Angle(arctan2((cos(epsilon) * cos(self.dec) * sin(self.ra) + sin(epsilon) * sin(self.dec)), (cos(self.dec) * cos(self.ra))), u.rad)

    @property
    def b_rate(self):
        num = self.dec_rate * (cos(epsilon) * cos(self.dec) + sin(epsilon) * sin(self.dec) * sin(self.ra)) - self.ra_rate * cos(self.dec) * cos(self.ra) * sin(epsilon)
        denom = sqrt(1 - (cos(epsilon) * sin(self.dec) - cos(self.dec) * sin(self.ra) * sin(epsilon))**2)
        return num / denom

    @property
    def l_rate(self):
        num = self.ra_rate * (1/cos(self.ra)) * (cos(epsilon) * (1/cos(self.ra)) + sin(epsilon) * tan(self.ra))
        denom = 1 + ((1/cos(self.ra)) * sin(epsilon) + cos(epsilon) * tan(self.ra))**2
        return num / denom




    def estimate_mag(self):
        '''
        Estimate the apparent magnitude of a TNO
        https://iopscience.iop.org/article/10.3847/1538-3881/ab18a9/pdf for light curves
        '''
        t = ts.tdb(jd=self.epoch.tdb.jd)
        e = earth.at(t)
        s = sun.at(t)
        sx, sy, sz = s.ecliptic_xyz().au
        ex, ey, ez = e.ecliptic_xyz().au
        earth_dist = ((ex-sx)**2 + (ey-sy)**2 + (ez-sz)**2)**0.5

        q = (self.r_helio.au**2 + self.delta.au**2 - earth_dist)/(2 * self.r_helio.au * self.delta.au)

        ## pyephem
        beta = np.arccos(q)
        beta[np.where(q <= -1)[0]] = np.pi * u.rad
        beta[np.where(q >= 1)[0]] = 0 * u.rad

        Psi_1 = np.exp(-3.332 * np.tan(beta/2)**0.631)
        Psi_2 = np.exp(-1.862 * np.tan(beta/2)**1.218)
        mag = self.H + 5 * np.log10(self.r_helio.au * self.delta.au)# / earth_dist**2)


        not_zero = np.where((Psi_1 != 0) | (Psi_2 != 0))[0]
        mag[not_zero] -= 2.5 * np.log10((1 - self.G[not_zero]) * Psi_1[not_zero] + self.G[not_zero] * Psi_2[not_zero])

        return mag



#    def plot_radec(self, color='black', alpha=0.5, zoom=False, galactic_plane=False, ecliptic_plane=True):
#        '''
#        Plot the right ascension and declination of each object on a Mollweide
#        projection. (See https://en.wikipedia.org/wiki/Mollweide_projection)
#        '''
#        if zoom == True:
#
#            fig = plt.figure(figsize=(12, 8))
#            ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
#
#            xdata = self.ra.degree
#            xdata[xdata > 180] -= 360
#
#            xmin=np.min(xdata)
#            xmax=np.max(xdata)
#            ymin=np.min(self.dec.degree)
#            ymax=np.max(self.dec.degree)
#
#            ax.scatter(-xdata, self.dec.degree, color='black', alpha=0.5)
#
#            xticks = np.linspace(-xmax, -xmin, 8)
#            yticks = np.linspace(ymin, ymax, 8)
#
#            ax.set_xticks(xticks)
#            ax.set_yticks(yticks)
#            xticklabels = [r'${:.2f}\degree$'.format(-value) for value in xticks]
#            yticklabels = [r'${:.2f}\degree$'.format(value) for value in yticks]
#            ax.set_xticklabels(xticklabels)
#            ax.set_yticklabels(yticklabels)
#
#            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                              linewidth=1, color='gray', alpha=0.5, linestyle='-')
#            gl.xlocator = mticker.FixedLocator(xticks)
#            gl.ylocator = mticker.FixedLocator(yticks)
#
#            gl.bottom_labels = False
#            gl.top_labels = False
#            gl.left_labels = False
#            gl.right_labels = False
#
#            xrange = xmax - xmin
#            yrange = ymax - ymin
#
#            try:
#                ax.set_extent([-xmax - xrange * 0.05, -xmin + xrange * 0.05,
#                               ymin - yrange * 0.05, ymax + yrange * 0.05], crs=ccrs.PlateCarree())
#            except:
#                ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
#
#            if ecliptic_plane == True:
#                def radec2project(ra, dec):
#                    ra[ra>180] -= 360
#                    return (ra, dec)
#
#                plane_lon = np.linspace(-np.pi, np.pi, 1000)
#                plane_lat = np.zeros(len(plane_lon))
#                ecl_plane = np.zeros([len(plane_lat), 2])
#
#                for i in range(len(plane_lat)):
#                    ecl_plane[i] = ephem.Equatorial(ephem.Ecliptic(plane_lon[i], plane_lat[i])).get()
#
#                x, y = radec2project(np.degrees(ecl_plane.T[0]), np.degrees(ecl_plane.T[1]))
#                ax.plot(x[1:999], -y[1:999], 'r-', zorder=2)
#
#            if galactic_plane == True:
#
#                galactic = SkyCoord(l=np.linspace(0, 360, 1000)*u.degree, b=np.zeros(1000)*u.degree, frame='galactic')
#                galxdata = galactic.icrs.ra.degree
#                galxdata[galxdata > 180] -= 360
#                galydata = galactic.icrs.dec.degree
#                order = galxdata.argsort()
#                ax.plot(-galxdata[order][1:999], -galydata[order][1:999], 'g-', zorder=3)
#
#        else:
#
#            fig = plt.figure(figsize=(12, 8))
#            ax = fig.add_subplot(111, projection='mollweide')
#            ax.grid(True)
#
#            xdata = self.ra
#            xdata[xdata.value > np.pi] -= 2*np.pi * u.rad
#            ax.set_xticklabels([r'$150\degree$', r'$120\degree$', r'$90\degree$' ,
#                                r'$60\degree$' , r'$30\degree$' , r'$0\degree$'  ,
#                                r'$330\degree$', r'$300\degree$', r'$270\degree$',
#                                r'$240\degree$', r'$210\degree$'])
#
#            ax.scatter(-xdata, self.dec, color=color, alpha=alpha)
#
#            if ecliptic_plane == True:
#                def radec2project(ra, dec):
#                    ra[ra>180] -= 360
#                    return (ra, dec)
#
#                plane_lon = np.linspace(-np.pi, np.pi, 1000)
#                plane_lat = np.zeros(len(plane_lon))
#                ecl_plane = np.zeros([len(plane_lat), 2])
#
#                for i in range(len(plane_lat)):
#                    ecl_plane[i] = ephem.Equatorial(ephem.Ecliptic(plane_lon[i], plane_lat[i])).get()
#
#                x, y = radec2project(np.degrees(ecl_plane.T[0]), np.degrees(ecl_plane.T[1]))
#                ax.plot(np.radians(x[1:999]), -np.radians(y[1:999]), 'r-', zorder=2)
#
#            if galactic_plane == True:
#
#                galactic = SkyCoord(l=np.linspace(0, 360, 1000)*u.degree, b=np.zeros(1000)*u.degree, frame='galactic')
#                galxdata = galactic.icrs.ra.rad
#                galxdata[galxdata > np.pi] -= 2*np.pi
#                galydata = galactic.icrs.dec.rad
#                order = galxdata.argsort()
#                ax.plot(-galxdata[order][1:999], galydata[order][1:999], 'g-', zorder=3)
#
#        return fig, ax
#
