import os
import sys

import numpy as np
from numpy import pi, arcsin, arctan2, log10, prod

import scipy
from scipy.interpolate import interpn

import pandas as pd
import healpy as hp # 2005ApJ...622..759G
import ephem
from numba import jit
from math import factorial, sqrt, sin, cos
import math

'''
All angle objects are in RADIANS!!!
'''

G = 0.15

# Change this to wherever you have the data
HPIX_DIR = '/Users/kjnapier/Research/TNOs/ClusteringOfETNOs/Data/hpix_data'

# Discovery pointings
DES_hpix = pd.read_csv(os.path.join(HPIX_DIR,'DES_HPIX64.csv'))
DES_S = pd.read_csv(os.path.join(HPIX_DIR,'DES_SN_S_HPIX1024.csv'))
DES_X12 = pd.read_csv(os.path.join(HPIX_DIR,'DES_SN_X12_HPIX1024.csv'))
DES_X3 = pd.read_csv(os.path.join(HPIX_DIR,'DES_SN_X3_HPIX1024.csv'))
DES_C12 = pd.read_csv(os.path.join(HPIX_DIR,'DES_SN_C12_HPIX1024.csv'))
DES_C3 = pd.read_csv(os.path.join(HPIX_DIR,'DES_SN_C3_HPIX1024.csv'))
DES_E = pd.read_csv(os.path.join(HPIX_DIR,'DES_SN_E_HPIX1024.csv'))
DES_SN_shallow = pd.concat([DES_S, DES_X12, DES_C12, DES_E])
DES_SN_deep = pd.concat([DES_X3, DES_C3])

# These are the HPIX that I fixed
OSSOS_hpix = pd.read_csv(os.path.join(HPIX_DIR,'Kev_OSSOS_HPIX.csv'))

B15BS = OSSOS_hpix[OSSOS_hpix.blockname == '15BS']
B15BT = OSSOS_hpix[OSSOS_hpix.blockname == '15BT']
B13BL = OSSOS_hpix[OSSOS_hpix.blockname == '13BL']
B14BH = OSSOS_hpix[OSSOS_hpix.blockname == '14BH']
B15BC = OSSOS_hpix[OSSOS_hpix.blockname == '15BC']
B15BC_2 = OSSOS_hpix[OSSOS_hpix.blockname == '15BC_2']
B15BD = OSSOS_hpix[OSSOS_hpix.blockname == '15BD']
B15AP = OSSOS_hpix[OSSOS_hpix.blockname == '15AP']
B13AE = OSSOS_hpix[OSSOS_hpix.blockname == '13AE']
B15AM = OSSOS_hpix[OSSOS_hpix.blockname == '15AM']
B13AO = OSSOS_hpix[OSSOS_hpix.blockname == '13AO']

ctio_hpix = pd.read_csv(os.path.join(HPIX_DIR,'ctio_hpix.csv'))
kpno_hpix = pd.read_csv(os.path.join(HPIX_DIR,'kpno_hpix.csv'))
magellan_hpix = pd.read_csv(os.path.join(HPIX_DIR,'magellan_hpix.csv'))
subaru_hpix = pd.read_csv(os.path.join(HPIX_DIR,'subaru_hpix.csv'))
ctio_goblin = pd.read_csv(os.path.join(HPIX_DIR,'ctio_goblin.csv'))
lbt_goblin = pd.read_csv(os.path.join(HPIX_DIR,'lbt_goblin.csv'))
subaru_goblin = pd.read_csv(os.path.join(HPIX_DIR,'subaru_goblin.csv'))
magellan_goblin = pd.read_csv(os.path.join(HPIX_DIR,'magellan_goblin.csv'))

ctio_st = pd.concat([ctio_hpix, ctio_goblin])
lbt_st = lbt_goblin
magellan_st = pd.concat([magellan_hpix, magellan_goblin])
subaru_st = pd.concat([subaru_hpix, subaru_goblin])
kpno_st = kpno_hpix


from skyfield.api import Topos, Loader
import ephem

# Specify single epoch for detections.
epoch = ephem.julian_date(ephem.date('2014/01/01'))

# Compute cartesian position of Earth at given epoch.
load = Loader('./Skyfield-Data', expire=False)
planets = load('de423.bsp')
earth = planets['earth']
ts = load.timescale()
t = ts.tai(jd=epoch+0.000428) #37 leap seconds
x_earth, y_earth, z_earth = earth.at(t).position.au # earth ICRS position


epsilon =  23.43929111 * pi/180. # obliquity, hard-coded
cosepsilon = cos(epsilon)
sinepsilon = sin(epsilon)


def check_if_in_survey(df, survey_hpix, maglim, NSIDE, ST=False):
    if NSIDE == 64:
        return df[df.hpix64.isin(survey_hpix.HPIX_64.values).values & 
                  (df.mag < maglim).values & (df.peri > 30).values]
    elif NSIDE == 1024:
        return df[df.hpix1024.isin(survey_hpix.HPIX_1024.values).values & 
                  (df.mag < maglim).values & (df.peri > 30).values]       
    elif ST == True:
        if NSIDE == 64:
            return df[df.hpix64.isin(survey_hpix.HPIX_64.values).values & 
                      (df.mag < maglim).values & (df.peri > 30).values &
                      (df.sun_distance > 50).values]
        elif NSIDE == 128:
            return df[df.hpix128.isin(survey_hpix.HPIX_128.values).values & 
                      (df.mag < maglim).values & (df.peri > 30).values &
                      (df.sun_distance > 50).values]
    

"""
Calculate eccentric anomaly using fixed point iteration 
(just like Kepler did in 1621)
"""
@jit
def calcE(e, ma, N):
    E = ma
    for kk in range(N):
        E = ma + e * np.sin(E)
    return E


@jit
# Calculate the RA and DEC on the epoch date. Very simplified version of Ed's code
def calcradec(e, E, a, omega, Omega, inc, H):
    cosarg = np.cos(omega)
    sinarg = np.sin(omega)
    cosnode = np.cos(Omega)
    sinnode = np.sin(Omega)
    cosi = np.cos(inc)
    sini = np.sin(inc)
    
    # compute true anomaly v
    v = 2 * arctan2((1 + e)**0.5 * np.sin(E/2), (1 - e)**0.5 * np.cos(E/2))
    # compute the distance to the central body r
    r = a * (1 - e * np.cos(E))
    
    # obtain the position o vector
    ox = r * np.cos(v)
    oy = r * np.sin(v)
        
    X0 = (ox * (cosarg * cosnode - sinarg * sinnode * cosi)
          - oy * (sinarg * cosnode + cosarg * sinnode * cosi))
    Y0 = (ox * (cosarg * sinnode + sinarg * cosnode * cosi)
          + oy * (cosarg * cosnode * cosi - sinarg * sinnode))
    Z0 = ox * (sinarg * sini) + oy * (cosarg * sini)
    r = (X0**2 + Y0**2 + Z0**2)**0.5
    
    # transfer ecliptic to ICRS and shift to Geocentric (topocentric)
    X = X0 - x_earth 
    Y = Y0 * cosepsilon - Z0 * sinepsilon - y_earth
    Z = Y0 * sinepsilon + Z0 * cosepsilon - z_earth
    delta = (X**2 + Y**2+ Z**2)**0.5
    # Cartesian to spherical coordinate
    dec = arcsin(Z/delta)
    ra = arctan2(Y, X) % (2 * pi)
    mag2 = H + 2.5 * log10(r**2 * delta**2)

    c = (r**2 + delta**2 - 1)/(2*r*delta)

    beta = np.zeros(len(c))
    beta[np.where(c <= 1)[0]] = np.pi
    beta[np.where(c >= 1)[0]] = 0

    tb2 = np.tan(beta/2.0);

    psi_t = tb2**0.63
    Psi_1 = np.exp(-3.33*psi_t)

    psi_t = tb2**1.22
    Psi_2 = np.exp(-1.87*psi_t)
    mag = H + 5.0*np.log10(r*delta)
    not_zero = np.where((Psi_1 != 0) | (Psi_2 != 0))[0]
    mag[not_zero] -= 2.5*np.log10((1-G)*Psi_1[not_zero] + G*Psi_2[not_zero])

    return ra, dec, delta, r, mag, mag2

# Convert RA, DEC to healpix index
def radec_to_index(ra, dec, NSIDE=64, nest=True):
    return hp.pixelfunc.ang2pix(NSIDE, -dec+pi/2, ra, nest)

# Calculate the likelihood of our detections in a given outer Solar System model
def calcL(N_iter, N_obj, varpiGrid, OmegaGrid,
          hist, varpiObj, OmegaObj, varpi, Omega):
    
    f = interp2d(varpiGrid, OmegaGrid, hist, kind='linear')

    joined = np.array([varpi, Omega]).T
    
    simprobs = []
    for kk in range(N_iter):
        prob = []
        samplepoints = np.random.randint(0, len(varpi), N_obj)
        samples = joined[samplepoints]
        for coordinate in samples:
            prob.append(f(coordinate[0], coordinate[1])[0])
        simprobs.append(prod(prob))

    ourprob = []
    for coord in zip(varpiObj, OmegaObj):
        ourprob.append(f(coord[0], coord[1])[0])
        
    ourprob = prod(ourprob)

    likelihood = len(np.where(simprobs < ourprob)[0])/len(simprobs)
    
    return likelihood, simprobs


# def mod_varpi(omega, Omega):
#     varpi = (omega + Omega) % 360
#     varpi = varpi - 360 * (varpi > 180)
#     return varpi

# def mod_Omega(in_Omega):
#     Omega = in_Omega - 360 * (in_Omega > 180)
#     return Omega

# Functions defining the orthonormal basis.
# All angular elements must be expressed in radians!
# x, y, p, and q accept and return numpy float arrays.

@jit
def Gamma(e):
    return 1 - np.array(list(map(sqrt, 1 - e**2)))

@jit
def Z(e, i):
    return np.array(list(map(sqrt, 1 - e**2))) * (1 - np.array(list(map(cos, i))))

@jit
def x(e, varpi):
    return np.array(list(map(sqrt, 2 * Gamma(e)))) * np.array(list(map(cos, varpi)))

@jit
def y(e, varpi):
    return np.array(list(map(sqrt, 2 * Gamma(e)))) * np.array(list(map(sin, varpi)))

@jit
def p(e, i, node):
    return np.array(list(map(sqrt, 2 * Z(e, i)))) * np.array(list(map(cos, node)))

@jit
def q(e, i, node):
    return np.array(list(map(sqrt, 2 * Z(e, i)))) * np.array(list(map(sin, node)))

def s_varpi(sv):
    return (sv.omega + sv.node) % (2*np.pi)
