#!/usr/bin/env python

from __future__ import division

import pyOrbfit as orbfit
import numpy as np
import ephem
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Ellipse

class Orbfit(object):
    """Class for interface to orbit-fitting code"""

    def __init__(self, *args, **kwargs):
        """
        Takes a set of input observations and fits the orbital elements. These can then be used to
        compute the position and error ellipse for the object on an arbitrary date.
        Requires: observatories.dat and binEphem.423.

        keyword arguments:
        obsfile -- file containing observation data, either in MPC format or as
                   JD HH:MM:SS.SS +DD:MM:SS.SS ERROR OBSCODE  or as
                   YYYY MM DD.DDDDD HH:MM:SS.SS +DD:MM:SS.SS ERROR OBSCODE
                   where the error is the astrometric error on each measurement in arcsec (default 0.2)
                   and OBSCODE is the observatory code in observatories.dat
                   The number of digits in each field is unrestricted, but ra and dec must not contain spaces.
        The following options should be used only if obsfile is None. Then they must all be used. If obsfile
        is supplied then these arguments are ignored.
        dates -- list of observation dates, either in JD or YYYY MM DD.DDDDD format
        ra -- list of RA in HH:MM:SS.SS format.
        dec -- list of DEC in +DD:MM:SS.SS format
        obscode -- list of observatory codes
        err -- measurement error in arcsec

        Methods:
           get_elements() -- returns aei orbital elements and their errors
           get_elements_abg() -- returns abg orbital elements and their errors
           barycentric_distance() -- barycentric distance and error
           perihelion() -- perihelion and error
           aphelion() -- aphelion and error
           cov_pq() -- covariance matrix for perihelion, aphelion
           plotEllipse() -- generates an error ellipse for a given covariance matrix
           predict_pos() -- predicted ra, dec, and error ellipse on a given date
           ellipticalBody() -- returns an EllipticalBody object with these orbital parameters, suitable for use by pyEphem
        """

        if 'file' in kwargs.keys():    # observations supplied in a file
            obsfile = kwargs['file']
            with open(obsfile, 'r') as fobs:
                self.nobs = 0 # determine the number of observations
                for line in fobs:
                    if line[0] != '#':
                        self.nobs+=1
                fobs.seek(0)   # rewind to the beginning
                self.obsarray = orbfit.OBSERVATION_ARRAY(self.nobs)   # create empty array of observations
                iobs = 0
                for line in fobs:
                    if line[0] != '#':
                        thisobs = orbfit.OBSERVATION()
                        orbfit.scan_observation(line, thisobs)
                        orbfit.add_to_obsarray(self.obsarray, iobs, thisobs)
                        iobs+=1
                if 'err' in kwargs.keys(): # change: obserr need not be constant
                    obserr = kwargs['err']
        elif len(args):              # observations supplied in a Catalog object
            self.nobs = len(args[0])
            self.obsarray = orbfit.OBSERVATION_ARRAY(self.nobs)
            for iobs, pt in enumerate(args[0]):
                thisobs = orbfit.OBSERVATION()
                if pt.err is None: pt.err=0.15           # This is to avoid silent crash.
                if pt.obscode is None: pt.obscode=807    # Ditto.
                obsline = str(ephem.julian_date(pt.date))+' '+str(pt.ra)+' '+str(pt.dec)+' '+str(pt.err)+' '+str(pt.obscode)
                orbfit.scan_observation(obsline, thisobs)
                orbfit.add_to_obsarray(self.obsarray, iobs, thisobs)
        else:                           # observations specified in input lists
            required_keys = ['dates', 'ra', 'dec', 'obscode']
            for k in required_keys:
                if k not in kwargs.keys():
                    raise KeyError('keyword '+k+' is missing')
            obsdate = kwargs['dates']
            ra = kwargs['ra']
            dec = kwargs['dec']
            obscode = kwargs['obscode']
            self.nobs = len(obsdate)
            if 'err' in kwargs.keys(): # change: obserr need not be constant
                obserr = kwargs['err']
                if not np.iterable(obserr): obserr = [obserr for i in range(self.nobs)]
            else:
                obserr = [0.15 for i in range(self.nobs)]       # obervation astrometric error defaults to 0.15"
            assert self.nobs==len(ra) and self.nobs==len(dec) and self.nobs==len(obscode) and self.nobs==len(obserr)
            self.obsarray = orbfit.OBSERVATION_ARRAY(self.nobs)   # create empty array of observations
            for iobs in range(self.nobs):   # fill the OBSERVATION_ARRAY
                thisobs = orbfit.OBSERVATION()
                obsline = str(ephem.julian_date(obsdate[iobs]))+' '+str(ra[iobs])+' '+str(dec[iobs])+' '+str(obserr[iobs])+' '+str(obscode[iobs])
                orbfit.scan_observation(obsline, thisobs)
                orbfit.add_to_obsarray(self.obsarray, iobs, thisobs)

        # At this point we have filled obsarray with the observations. Now fit the orbit.

#        orbfit.set_ephem_file('/Users/gerdes/TNO/pyOrbfit/binEphem.430')      # pick up the correct ephemeris and observatories file.
#        orbfit.set_observatory_file('/Users/gerdes/TNO/pyOrbfit/observatories.dat')   # Will need to handle pathnames more elegantly.

        self.orbit_abg = orbfit.PBASIS()
        self.orbit_xyz = orbfit.XVBASIS()
        self.orbit_aei = orbfit.ORBIT()
        self.cov_abg = orbfit.dmatrix(1,6,1,6)   # need to make this one globally visible since it's needed by predict_pos()
        cov_xyz = orbfit.dmatrix(1,6,1,6)
        cov_aei = orbfit.dmatrix(1,6,1,6)
        derivs  = orbfit.dmatrix(1,6,1,6)
    # fittype is 6 if normal fit, 5 if energy constraint was used.
        self.fittype, self.chisq, self.ndof = orbfit.fit_observations(self.obsarray, self.nobs, self.orbit_abg, self.cov_abg, None)  # abg orbit elements
    #   Transform the orbit basis and get the deriv. matrix
        orbfit.pbasis_to_bary(self.orbit_abg, self.orbit_xyz, derivs)
        orbfit.orbitElements(self.orbit_xyz, self.orbit_aei) # aei orbit elements
#        self.orbit_xyz.jd0 = orbfit.cvar.jd0     # time zeropoint
# covariance matrices
        orbfit.covar_map(self.cov_abg, derivs, cov_xyz, 6, 6)    # map the covariance matrix to xyz basis
        orbfit.aei_derivs(self.orbit_xyz, derivs)   # Get partial derivative matrix from xyz to aei
        orbfit.covar_map(cov_xyz, derivs, cov_aei, 6, 6) # map covariance matrix from xyz to aei
# This is a hack to create matrices python can actually use. We have trouble wrapping double pointers.
        c = orbfit.doubleArray(36)
        orbfit.flatten_cov(self.cov_abg, 6, c)
        self.covar_abg = np.array([c[i] for i in range(36)]).reshape(6,6)   # a bona-fide numpy array.
        c = orbfit.doubleArray(36)  # It's necessary to reallocate the space in memory
        orbfit.flatten_cov(cov_xyz, 6, c)
        self.covar_xyz = np.array([c[i] for i in range(36)]).reshape(6,6)
        c = orbfit.doubleArray(36)
        orbfit.flatten_cov(cov_aei, 6, c)
        self.covar_aei = np.array([c[i] for i in range(36)]).reshape(6,6)
        self.elements, self.elements_errs = self.get_elements()
        self.jd0 = orbfit.cvar.jd0
        self.lat0 = orbfit.cvar.lat0
        self.lon0 = orbfit.cvar.lon0
        self.xBary = orbfit.cvar.xBary
        self.yBary = orbfit.cvar.yBary
        self.zBary = orbfit.cvar.zBary

    def mean_anomaly(self):
        GM = 4.*np.pi*np.pi/1.0000378 # solar gravitation, from orbfit
        combinedMass = GM * 1.00134  #Alter GM to account for total SS mass, from orbfit
        DAY = 1/365.25
        if np.isnan(self.elements['top']) or self.elements['a']<0:
            mu = -999
        else:
            mu = (self.jd0-self.elements['top'])*np.sqrt(combinedMass/self.elements['a']**3)*DAY*180/np.pi
            if mu<0: mu+=360
        return mu

    def barycentric_distance(self):
        """Barycentric distance and error"""
        xbary = orbfit.cvar.xBary
        ybary = orbfit.cvar.yBary
        zbary = orbfit.cvar.zBary
        dbary = np.sqrt(xbary**2 + ybary**2 + (zbary-1/self.orbit_abg.g)**2)
        dbary_err = dbary**2*np.sqrt(self.covar_abg[4][4])
        return dbary, dbary_err

    def perihelion(self):
        """Perihelion and error"""
        a = self.orbit_aei.a
        e = self.orbit_aei.e
        p = (1-e)*a
        p_err = np.sqrt(a**2*self.covar_aei[1][1] + (1-e)**2*self.covar_aei[0][0] + 2*a*(1-e)*self.covar_aei[0][1])
        return p, p_err

    def aphelion(self):
        """Aphelion and error"""
        a = self.orbit_aei.a
        e = self.orbit_aei.e
        q = (1+e)*a
        q_err = np.sqrt(a**2*self.covar_aei[1][1]+(1+e)**2*self.covar_aei[0][0] + 2*a*(1+e)*self.covar_aei[0][1])
        return q, q_err

    def cov_pq(self):
        """Covariance matrix for perihelion and aphelion (p,q)"""
        p, dp = self.perihelion()
        q, dq = self.aphelion()
        a = self.orbit_aei.a
        e = self.orbit_aei.e
        J = np.array([[1-e, -a],[1+e, a]])   # Jacobian of transformation
#       J = np.array([[(1-e)/a, -a/e],[(1+e)/a, a/e]])   # Jacobian of transformation
        cov_ae = np.array([[self.covar_aei[0][0], self.covar_aei[0][1]], [self.covar_aei[1][0], self.covar_aei[1][1]]])
        cpq = J.dot(cov_ae.dot(J.T))
        return cpq

    def get_elements(self):
        """Returns dictionaries of the standard aei orbital elements and their errors"""
        elements = {}
        elements_errs = {}
        elements['a'] = self.orbit_aei.a
        try:                                # semimajor axis (AU)
            elements_errs['a'] = np.sqrt(self.covar_aei[0][0])
        except RunTimeWarning:
            elements_errs['a'] = -99
        elements['e'] = self.orbit_aei.e
        try:                                # eccentricity
            elements_errs['e'] = np.sqrt(self.covar_aei[1][1])
        except RunTimeWarning:
            elements_errs['e'] = -99
        elements['i'] = self.orbit_aei.i
        try:                                # inclination (deg)
            elements_errs['i'] = np.sqrt(self.covar_aei[2][2])/(np.pi/180)
        except RunTimeWarning:
            elements_errs['i'] = -99
        elements['lan'] = self.orbit_aei.lan
        try:                             # longitude of ascending node (deg)
            elements_errs['lan'] = np.sqrt(self.covar_aei[3][3])/(np.pi/180)
        except RunTimeWarning:
            elements_errs['lan'] = -99
        elements['aop'] = self.orbit_aei.aop
        try:                            # argument of perihelion (deg)
            elements_errs['aop'] = np.sqrt(self.covar_aei[4][4])/(np.pi/180)
        except RunTimeWarning:
            elements_errs['aop'] = -99
        elements['top'] = self.orbit_aei.T
        try:                                # time of periapsis (JD)
            elements_errs['top'] = np.sqrt(self.covar_aei[5][5])/orbfit.DAY
        except RunTimeWarning:
            elements_errs['top'] = -99
        return elements, elements_errs

    def get_elements_abg(self):
        """Returns dictionaries of the abg orbital elements and their errors"""
        elements = {}
        elements_errs = {}
        elements['a'] = self.orbit_abg.a                                    # alpha, angular position at t=0
        elements_errs['a'] = np.sqrt(self.covar_abg[0][0])
        elements['adot'] = self.orbit_abg.adot                              # (dx/dt)/z at t=0
        elements_errs['adot'] = np.sqrt(self.covar_abg[1][1])
        elements['b'] = self.orbit_abg.b                                    # same as a, adot in y-direction
        elements_errs['b'] = np.sqrt(self.covar_abg[2][2])
        elements['bdot'] = self.orbit_abg.bdot
        elements_errs['bdot'] = np.sqrt(self.covar_abg[3][3])
        elements['g'] = self.orbit_abg.g                                    # 1/z of KBO at t=0
        elements_errs['g'] = np.sqrt(self.covar_abg[4][4])
        elements['gdot'] = self.orbit_abg.gdot                              # (dz/dt)/z at t = 0
        elements_errs['gdot'] = np.sqrt(self.covar_abg[5][5])
        return elements, elements_errs

    def plotEllipse(self, pos, cov, ind_x, ind_y, edge='k', face='b', alpha=0.5, scale_x=1.0, scale_y=1.0, label=None):
        """Returns an error ellipse, suitable for plotting, given a covariance matrix and the indices of the two variables
        Usage: ellipsePlot=plotEllipse([x,y], covar_aei, 0, 1, edge='b', face='b', alpha=0.7, label='Error ellipse')
        """
        cov2x2 = np.array([cov[ind_x][ind_x]/scale_x**2, cov[ind_x][ind_y]/(scale_x*scale_y), \
                           cov[ind_y][ind_x]/(scale_x*scale_y), cov[ind_y][ind_y]/scale_y**2]).reshape(2,2)
        U, s , Vh = np.linalg.svd(cov2x2)
        orient = np.arctan2(U[1,0],U[0,0])*180/np.pi
        ellipsePlot = Ellipse(xy=pos, width=2.0*np.sqrt(s[0]), height=2.0*np.sqrt(s[1]), angle=orient, facecolor=face, edgecolor=edge, alpha=alpha, label=label)
        return ellipsePlot

    def predict_pos(self, date, obscode=807):
        """
        Computes ra, dec, and error ellipse at the date(s) and observatory specified.
        Date is an ephem.date object.
        Returns a corresponding dictionary with keywords
            'ra' -- ra in ICRS coords, returned as an ephem.Angle object
            'dec' -- dec in ICRS coords, ditto
            'err' -- error ellipse (see below)
            'elong' -- solar elongation in degrees
            'opp' -- opposition angle in degrees
        where the error ellipse, indexed by 'err', is itself a dictionary with keywords:
            'a'  -- semimajor axis in arcsec
            'b'  -- semiminor axis in arcsec
            'PA' -- position angle in degrees, north through east

        """
        p_in = self.orbit_abg
        jd0 = orbfit.cvar.jd0     # zeropoint of time scale
        # Create space for various useful matrices
        sigxy = orbfit.dmatrix(1,2,1,2)
        derivs = orbfit.dmatrix(1,2,1,2)
        covecl = orbfit.dmatrix(1,2,1,2)
        coveq = orbfit.dmatrix(1,2,1,2)

        # Fill the OBSERVATION structure
        futobs = orbfit.OBSERVATION()
        futobs.obscode = obscode
        futobs.obstime = (ephem.julian_date(date)-jd0)*orbfit.DAY
        futobs.xe = -999  # force evaluation of earth3D
        dx = orbfit.dvector(1,6)
        dy = orbfit.dvector(1,6)
        # Sometimes hangs in next line... why?
        thetax, thetay = orbfit.kbo2d(p_in, futobs, dx, dy)
        # Predicted position, in abg basis:
        orbfit.predict_posn(p_in, self.cov_abg, futobs, sigxy)
        solar_elongation = orbfit.elongation(futobs)/orbfit.DTOR       # solar elongation in degrees
        opposition_angle = orbfit.opposition_angle(futobs)/orbfit.DTOR # opposition angle in degrees
        lat_ec, lon_ec = orbfit.proj_to_ec(futobs.thetax, futobs.thetay, orbfit.cvar.lat0, orbfit.cvar.lon0, derivs)  # project to ecliptic coords
        orbfit.covar_map(sigxy, derivs, covecl, 2, 2)    # map the covariance
        ra_eq, dec_eq = orbfit.ec_to_eq(lat_ec, lon_ec, derivs)    # transform to ICRS to compute ra, dec
        if ra_eq<0: ra_eq += 2*np.pi        # convert angle to (0,2pi) range
        orbfit.covar_map(covecl, derivs, coveq, 2, 2)    # map the covariance

        # Compute ICRS error ellipse
        c = orbfit.doubleArray(4)   # stoopid workaround for double pointers...
        orbfit.flatten_cov(coveq, 2, c)
        covar_eq = np.array([c[i] for i in range(4)]).reshape(2,2)
        xx = covar_eq[0][0]*np.cos(dec_eq)**2
        xy = covar_eq[0][1]*np.cos(dec_eq)
        yy = covar_eq[1][1]
        pos_angle = 0.5*np.arctan2(2.*xy,(xx-yy)) * 180./np.pi
        pos_angle = 90 - pos_angle     # convert to astronomy convention of measuring position angle North through East
        bovasqrd  = (xx + yy - np.sqrt((xx-yy)**2 + (2*xy)**2)) / (xx + yy + np.sqrt((xx-yy)**2 + (2*xy)**2))
        det = xx*yy-xy*xy
        a = (det/bovasqrd)**(1/4)/orbfit.ARCSEC   # semimajor, minor axes of error ellipse, in arcsec
        b = (det*bovasqrd)**(1/4)/orbfit.ARCSEC
        err_ellipse = dict(a=a, b=b, PA=pos_angle)   # store as a dictionary

        pos = dict(ra=ephem.hours(ra_eq), dec=ephem.degrees(dec_eq), err=err_ellipse, elong=solar_elongation, opp=opposition_angle)
        return pos

    def ellipticalBody(self, name=None):
        '''
        Packs the orbital parameters into an EllipticalBody object suitable for use by pyEphem
        '''
        DJD = 2415020   # Dublin Julian Date t=0 (used by pyephem Date object)
        b = ephem.EllipticalBody()
        if name is not None:
            b.name=name
        b._a = self.orbit_aei.a
        b._inc = self.orbit_aei.i
        b._Om = self.orbit_aei.lan
        b._om = self.orbit_aei.aop
        b._e = self.orbit_aei.e
        b._epoch = ephem.date('2000/01/01')    # J2000
        b._M = self.mean_anomaly()     # mean anomaly from perihelion (deg)
        b._epoch_M = ephem.date(orbfit.cvar.jd0-DJD)   # date for measurement of mean anomaly _M
        return b

    def orbit2df(self):
        '''
        Insert the orbit details into a pandas dataframe.
        '''
        orbit_cols = ['chisq', 'ndof', 'a', 'e', 'inc',
        'aop', 'node', 'peri_jd', 'peri_date', 'epoch_jd',
        'mean_anomaly', 'period', 'period_err',
        'a_err', 'e_err', 'inc_err', 'aop_err', 'node_err', 'peri_err', 'lat0', 'lon0',
        'xbary','ybary','zbary',
        'abg_a', 'abg_b', 'abg_g',
        'abg_adot', 'abg_bdot', 'abg_gdot',
        'abg_a_err', 'abg_b_err', 'abg_g_err',
        'abg_adot_err', 'abg_bdot_err', 'abg_gdot_err']
        index = np.arange(1)
        df = pd.DataFrame(columns=orbit_cols, index=index)

        elements, errs = self.get_elements()
        elements_abg, errs_abg = self.get_elements_abg()
        epoch = self.jd0

        peri_jd = -999 if np.isnan(elements['top']) else elements['top']
        peri_date = -999 if np.isnan(elements['top']) else str(ephem.date(elements['top']-2415020))
        peri_jd_err = -999 if (np.isnan(elements['top']) or np.isnan(errs['top'])) else errs['top']
        period = elements['a']**1.5 if elements['a']>0 else -999
        period_err = 3*np.sqrt(elements['a'])/2*errs['a'] if elements['a']>0 else -999

        orbit_data = np.array([round(self.chisq,2), self.ndof, round(elements['a'],4), round(elements['e'],6), round(elements['i'],6),
            round(elements['aop'],2), round(elements['lan'],3), round(peri_jd,2), peri_date, round(epoch,2),
            self.mean_anomaly(), round(period,3) , round(period_err,3),
            round(errs['a'],4), round(errs['e'],6), round(errs['i'],6), round(errs['aop'],2), round(errs['lan'],3), round(peri_jd_err,2), round(self.lat0,4), round(self.lon0,4),
            round(self.xBary,4), round(self.yBary,4), round(self.zBary,4),
            elements_abg['a'], elements_abg['b'], elements_abg['g'],
            elements_abg['adot'], elements_abg['bdot'], elements_abg['gdot'],
            errs_abg['a'], errs_abg['b'], errs_abg['g'],
            errs_abg['adot'], errs_abg['bdot'], errs_abg['gdot']
            ])
        df.ix[0] = orbit_data
        return df

    def residuals(self, obsdf, date_col='date', ra_col='ra', dec_col='dec'):
        '''
        Returns the residuals (in arcseconds) for the observations in obscat as a numpy array. Ordering is same as time-order of the observations,
        from earliest to latest.
        obsdf is a pandas dataframe containing the observations.
        date, ra, dec are in the column names specified.
        date is an ephem.date object
        ra, dec are ephem.Angle objects (typically hours, degrees)
        '''
        resids = []
        obsdf.sort_values('expnum', ascending=True, inplace=True)
        for ind, obs in obsdf.iterrows():
            if obs[date_col][0] != '#':
                pos_pred = self.predict_pos(ephem.date(obs[date_col])+0.53*ephem.second + float(obs['exptime'])*ephem.second/2)
                ra_pred, dec_pred = pos_pred['ra'], pos_pred['dec']
                sep = ephem.separation( (ephem.hours(obs[ra_col]),ephem.degrees(obs[dec_col])), (ra_pred, dec_pred))
                resid = sep*(180/np.pi)*3600    # convert to arcseconds
                resids.append(resid)
        return np.array(resids)




def main():
    pass

if __name__ == '__main__':
    main()
