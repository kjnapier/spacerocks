import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np

from astropy import units as u
from astropy.coordinates import Angle, Distance

from . import clibspacerocks

clibspacerocks.py_kepM_to_xyz.argtypes = [ctypes.c_int,
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_kepM_to_xyz.restype = ctypes.POINTER(ctypes.c_double)

def kepM_to_xyz(a, e, inc, arg, node, M):

    N = len(a)
    rock = clibspacerocks.py_kepM_to_xyz(N, a, e, inc, arg, node, M)
    arr = np.ctypeslib.as_array(rock, (6 * N,))

    x, y, z, vx, vy, vz = arr.reshape(N, 6).T

    x = Distance(x, u.au, allow_negative=True)
    y = Distance(y, u.au, allow_negative=True)
    z = Distance(z, u.au, allow_negative=True)
    
    return x, y, z, vx * u.au/u.d, vy * u.au/u.d, vz * u.au/u.d

clibspacerocks.py_kepE_to_xyz.argtypes = [ctypes.c_int,
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_kepE_to_xyz.restype = ctypes.POINTER(ctypes.c_double)

def kepE_to_xyz(a, e, inc, arg, node, E):

    N = len(a)
    rock = clibspacerocks.py_kepE_to_xyz(N, a, e, inc, arg, node, E)
    arr = np.ctypeslib.as_array(rock, (6 * N,))

    x, y, z, vx, vy, vz = arr.reshape(N, 6).T

    x = Distance(x, u.au, allow_negative=True)
    y = Distance(y, u.au, allow_negative=True)
    z = Distance(z, u.au, allow_negative=True)
    
    return x, y, z, vx * u.au/u.d, vy * u.au/u.d, vz * u.au/u.d


clibspacerocks.py_calc_kep_from_xyz.argtypes = [ctypes.c_int,
                                                ctypes.c_double,
                                                ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_calc_kep_from_xyz.restype = ctypes.POINTER(ctypes.c_double)

def calc_kep_from_xyz(mu, x, y, z, vx, vy, vz):

    N = len(x)
    rock = clibspacerocks.py_calc_kep_from_xyz(N, mu, x, y, z, vx, vy, vz)
    arr = np.ctypeslib.as_array(rock, (6 * N,))

    a, e, inc, arg, node, f = arr.reshape(N, 6).T
    a = Distance(a, u.au, allow_negative=True)
    e = np.ascontiguousarray(e)
    inc = Angle(inc, u.rad)
    arg = Angle(arg, u.rad)
    node = Angle(node, u.rad)
    f = Angle(f, u.rad)

    return a, e, inc, arg, node, f


clibspacerocks.py_calc_vovec_from_kep.argtypes = [ctypes.c_int,
                                                  ctypes.c_double,
                                                  ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                  ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                  ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                  ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_calc_vovec_from_kep.restype = ctypes.POINTER(ctypes.c_double)


def calc_vovec_from_kep(mu, a, e, r, E):

    N = len(a)
    rock = clibspacerocks.py_calc_vovec_from_kep(N, mu, a, e, r, E)
    arr = np.ctypeslib.as_array(rock, (3 * N,))

    vx, vy, vz = arr.reshape(N, 3).T
    return vx * u.au/u.d, vy * u.au/u.d, vz * u.au/u.d


clibspacerocks.py_calc_E_from_M.argtypes = [ctypes.c_int,
                                            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_calc_E_from_M.restype = ctypes.POINTER(ctypes.c_double)

def calc_E_from_M(e, M):

    N = len(e)
    rock = clibspacerocks.py_calc_E_from_M(N, e, M)
    E = np.ctypeslib.as_array(rock, (N,))

    return Angle(E, u.rad)


clibspacerocks.py_calc_M_from_E.argtypes = [ctypes.c_int,
                                            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_calc_M_from_E.restype = ctypes.POINTER(ctypes.c_double)

def calc_M_from_E(e, E):

    N = len(e)
    rock = clibspacerocks.py_calc_M_from_E(N, e, E)
    M = np.ctypeslib.as_array(rock, (N,))

    return Angle(M, u.rad)


clibspacerocks.py_calc_E_from_f.argtypes = [ctypes.c_int,
                                                       ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                       ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_calc_E_from_f.restype = ctypes.POINTER(ctypes.c_double)

def calc_E_from_f(e, f):

    N = len(e)
    rock = clibspacerocks.py_calc_E_from_f(N, e, f)
    E = np.ctypeslib.as_array(rock, (N,))

    return Angle(E, u.rad)


clibspacerocks.py_calc_f_from_E.argtypes = [ctypes.c_int,
                                                       ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                       ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_calc_f_from_E.restype = ctypes.POINTER(ctypes.c_double)

def calc_f_from_E(e, E):

    N = len(e)
    rock = clibspacerocks.py_calc_f_from_E(N, e, E)
    f = np.ctypeslib.as_array(rock, (N,))

    return Angle(f, u.rad)
