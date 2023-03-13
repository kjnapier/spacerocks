import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np

from astropy import units as u
from astropy.coordinates import Angle, Distance

from . import clibspacerocks

M_TO_AU = (u.m).to(u.au)
M_TO_KM = (u.m).to(u.km)

clibspacerocks.py_compute_topocentric_correction.argtypes = [ctypes.c_int,
                                                             ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                             ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                             ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                             ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_compute_topocentric_correction.restype = ctypes.POINTER(ctypes.c_double)

def compute_topocentric_correction(lat, lon, elevation, epoch):

    N = len(epoch)

    lat = np.atleast_1d(lat).astype(np.float64)
    lon = np.atleast_1d(lon).astype(np.float64)
    elevation = np.atleast_1d(elevation).astype(np.float64)

    if len(lat) < N:
        lat = np.repeat(lat, N)
        lon = np.repeat(lon, N)
        elevation = np.repeat(elevation, N)

    rock = clibspacerocks.py_compute_topocentric_correction(N, lat, lon, elevation, epoch)
    arr = np.ctypeslib.as_array(rock, (3 * N,)).copy()
    clibspacerocks.free_memory(rock)
    
    x, y, z = arr.reshape(N, 3).T
    
    return x * M_TO_KM, y * M_TO_KM, z * M_TO_KM

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
    arr = np.ctypeslib.as_array(rock, (6 * N,)).copy()
    clibspacerocks.free_memory(rock)
    
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
    arr = np.ctypeslib.as_array(rock, (6 * N,)).copy()
    clibspacerocks.free_memory(rock)
    
    x, y, z, vx, vy, vz = arr.reshape(N, 6).T

    x = Distance(x, u.au, allow_negative=True)
    y = Distance(y, u.au, allow_negative=True)
    z = Distance(z, u.au, allow_negative=True)
  
    
    return x, y, z, vx * u.au/u.d, vy * u.au/u.d, vz * u.au/u.d

clibspacerocks.py_kepf_to_xyz.argtypes = [ctypes.c_int,
                                          ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                          ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                          ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                          ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                          ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                          ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_kepf_to_xyz.restype = ctypes.POINTER(ctypes.c_double)

def kepf_to_xyz(a, e, inc, arg, node, f):

    N = len(a)
    rock = clibspacerocks.py_kepf_to_xyz(N, a, e, inc, arg, node, f)
    arr = np.ctypeslib.as_array(rock, (6 * N,)).copy()
    clibspacerocks.free_memory(rock)
    
    x, y, z, vx, vy, vz = arr.reshape(N, 6).T

    x = Distance(x, u.au, allow_negative=True)
    y = Distance(y, u.au, allow_negative=True)
    z = Distance(z, u.au, allow_negative=True)
  
    
    return x, y, z, vx * u.au/u.d, vy * u.au/u.d, vz * u.au/u.d


# clibspacerocks.py_calc_kep_from_xyz.argtypes = [ctypes.c_int,
#                                                 ctypes.c_double,
#                                                 ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
#                                                 ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
#                                                 ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
#                                                 ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
#                                                 ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
#                                                 ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

# clibspacerocks.py_calc_kep_from_xyz.restype = ctypes.POINTER(ctypes.c_double)

clibspacerocks.py_calc_kep_from_xyz.argtypes = [ctypes.c_int,
                                                ctypes.c_double,
                                                ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                                ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_calc_kep_from_xyz.restype = None

# def calc_kep_from_xyz(mu, x, y, z, vx, vy, vz):

#     N = len(x)
#     rock = clibspacerocks.py_calc_kep_from_xyz(N, mu, x, y, z, vx, vy, vz)
#     arr = np.ctypeslib.as_array(rock, (6 * N,)).copy()
#     clibspacerocks.free_memory(rock)
    
#     a, e, inc, arg, node, f = arr.reshape(N, 6).T
#     a = Distance(a, u.au, allow_negative=True)
#     e = np.ascontiguousarray(e)
#     inc = Angle(inc, u.rad)
#     arg = Angle(arg, u.rad)
#     node = Angle(node, u.rad)
#     f = Angle(f, u.rad)

#     return a, e, inc, arg, node, f

def calc_kep_from_xyz(mu, x, y, z, vx, vy, vz):

    N = len(x)
    arr = np.empty((6 * N,), dtype=np.float64)
    clibspacerocks.py_calc_kep_from_xyz(N, mu, x, y, z, vx, vy, vz, arr)
    
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
    rock = clibspacerocks.py_calc_vovec_from_kep(N, mu, a, e.astype(np.float64), r, E)
    arr = np.ctypeslib.as_array(rock, (3 * N,)).copy()
    clibspacerocks.free_memory(rock)

    vx, vy, vz = arr.reshape(N, 3).T
    
    return vx * u.au/u.d, vy * u.au/u.d, vz * u.au/u.d


clibspacerocks.py_calc_E_from_M.argtypes = [ctypes.c_int,
                                            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_calc_E_from_M.restype = ctypes.POINTER(ctypes.c_double)

def calc_E_from_M(e, M):

    N = len(e)
    rock = clibspacerocks.py_calc_E_from_M(N, e.astype(np.float64), M.astype(np.float64))
    E = np.ctypeslib.as_array(rock, (N,)).copy()
    clibspacerocks.free_memory(rock)

    return Angle(E, u.rad)


# clibspacerocks.py_calc_M_from_E.argtypes = [ctypes.c_int,
#                                             ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
#                                             ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

# clibspacerocks.py_calc_M_from_E.restype = ctypes.POINTER(ctypes.c_double)

# def calc_M_from_E(e, E):

#     N = len(e)
#     rock = clibspacerocks.py_calc_M_from_E(N, e, E)
#     # ptr = ctypes.cast(rock, ctypes.POINTER(ctypes.c_double * N))
#     # M = np.asarray(ptr.contents).copy() #
#     M = np.ctypeslib.as_array(rock, (N,)).copy()
#     clibspacerocks.free_memory(rock)

#     return Angle(M, u.rad)

clibspacerocks.py_calc_M_from_E.argtypes = [ctypes.c_int,
                                            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

#clibspacerocks.py_calc_M_from_E.restype = ctypes.POINTER(ctypes.c_double)

def calc_M_from_E(e, E):

    N = len(e)

    M = np.empty(N)
    clibspacerocks.py_calc_M_from_E(N, e, E, M)

    return Angle(M, u.rad)


clibspacerocks.py_calc_E_from_f.argtypes = [ctypes.c_int,
                                            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_calc_E_from_f.restype = ctypes.POINTER(ctypes.c_double)

def calc_E_from_f(e, f):

    N = len(e)
    rock = clibspacerocks.py_calc_E_from_f(N, e.astype(np.float64), f)
    E = np.ctypeslib.as_array(rock, (N,)).copy()
    clibspacerocks.free_memory(rock)
    return Angle(E, u.rad)


clibspacerocks.py_calc_f_from_E.argtypes = [ctypes.c_int,
                                            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_calc_f_from_E.restype = ctypes.POINTER(ctypes.c_double)

def calc_f_from_E(e, E):

    N = len(e)
    rock = clibspacerocks.py_calc_f_from_E(N, e.astype(np.float64), E.astype(np.float64))
    f = np.ctypeslib.as_array(rock, (N,)).copy()
    clibspacerocks.free_memory(rock)
    return Angle(f, u.rad)


clibspacerocks.py_correct_for_ltt.argtypes = [ctypes.c_int,
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]
clibspacerocks.py_correct_for_ltt.restype = ctypes.POINTER(ctypes.c_double)

clibspacerocks.py_correct_for_ltt_single_observer.argtypes = [ctypes.c_int,
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]
clibspacerocks.py_correct_for_ltt_single_observer.restype = ctypes.POINTER(ctypes.c_double)

# def correct_for_ltt(rocks, observers):

#     N = len(rocks)
#     a = rocks.a.au.astype(np.double)
#     e = rocks.e.astype(np.double)
#     inc = rocks.inc.rad.astype(np.double)
#     arg = rocks.arg.rad.astype(np.double)
#     node = rocks.node.rad.astype(np.double)
#     M = rocks.M.rad.astype(np.double)

#     ox = observers.x.au.astype(np.double)
#     oy = observers.y.au.astype(np.double)
#     oz = observers.z.au.astype(np.double)
#     ovx = observers.vx.to(u.au/u.day).value.astype(np.double)
#     ovy = observers.vy.to(u.au/u.day).value.astype(np.double)
#     ovz = observers.vz.to(u.au/u.day).value.astype(np.double)

#     rock = clibspacerocks.py_correct_for_ltt(N, a, e, inc, arg, node, M, ox, oy, oz, ovx, ovy, ovz)
#     arr = np.ctypeslib.as_array(rock, (6 * N,)).copy()
#     clibspacerocks.free_memory(rock)
  
#     x, y, z, vx, vy, vz = arr.reshape(N, 6).T

#     x = Distance(x, u.au, allow_negative=True)
#     y = Distance(y, u.au, allow_negative=True)
#     z = Distance(z, u.au, allow_negative=True)

#     return x, y, z, vx * u.au/u.day, vy * u.au/u.day, vz * u.au/u.day
    
# def correct_for_ltt_single_observer(rocks, observers):

#     N = len(rocks)
#     a = rocks.a.au.astype(np.double)
#     e = rocks.e.astype(np.double)
#     inc = rocks.inc.rad.astype(np.double)
#     arg = rocks.arg.rad.astype(np.double)
#     node = rocks.node.rad.astype(np.double)
#     M = rocks.M.rad.astype(np.double)

#     ox = observers.x.au.astype(np.double)
#     oy = observers.y.au.astype(np.double)
#     oz = observers.z.au.astype(np.double)
#     ovx = observers.vx.to(u.au/u.day).value.astype(np.double)
#     ovy = observers.vy.to(u.au/u.day).value.astype(np.double)
#     ovz = observers.vz.to(u.au/u.day).value.astype(np.double)

#     rock = clibspacerocks.py_correct_for_ltt(N, a, e, inc, arg, node, M, ox, oy, oz, ovx, ovy, ovz)
#     arr = np.ctypeslib.as_array(rock, (6 * N,)).copy()
#     clibspacerocks.free_memory(rock)
  
#     x, y, z, vx, vy, vz = arr.reshape(N, 6).T

#     x = Distance(x, u.au, allow_negative=True)
#     y = Distance(y, u.au, allow_negative=True)
#     z = Distance(z, u.au, allow_negative=True)

#     return x, y, z, vx * u.au/u.day, vy * u.au/u.day, vz * u.au/u.day

def correct_for_ltt(rocks, observers):

    N = len(rocks)
    x = rocks.x.au.astype(np.double)
    y = rocks.y.au.astype(np.double)
    z = rocks.z.au.astype(np.double)
    vx = rocks.vx.to(u.au/u.day).value.astype(np.double)
    vy = rocks.vy.to(u.au/u.day).value.astype(np.double)
    vz = rocks.vz.to(u.au/u.day).value.astype(np.double)

    ox = observers.x.au.astype(np.double)
    oy = observers.y.au.astype(np.double)
    oz = observers.z.au.astype(np.double)
    ovx = observers.vx.to(u.au/u.day).value.astype(np.double)
    ovy = observers.vy.to(u.au/u.day).value.astype(np.double)
    ovz = observers.vz.to(u.au/u.day).value.astype(np.double)

    rock = clibspacerocks.py_correct_for_ltt(N, x, y, z, vx, vy, vz, ox, oy, oz, ovx, ovy, ovz)
    arr = np.ctypeslib.as_array(rock, (6 * N,)).copy()
    clibspacerocks.free_memory(rock)
  
    x, y, z, vx, vy, vz = arr.reshape(N, 6).T

    x = Distance(x, u.au, allow_negative=True)
    y = Distance(y, u.au, allow_negative=True)
    z = Distance(z, u.au, allow_negative=True)

    return x, y, z, vx * u.au/u.day, vy * u.au/u.day, vz * u.au/u.day
    
def correct_for_ltt_single_observer(rocks, observers):

    N = len(rocks)
    x = rocks.x.au.astype(np.double)
    y = rocks.y.au.astype(np.double)
    z = rocks.z.au.astype(np.double)
    vx = rocks.vx.to(u.au/u.day).value.astype(np.double)
    vy = rocks.vy.to(u.au/u.day).value.astype(np.double)
    vz = rocks.vz.to(u.au/u.day).value.astype(np.double)

    ox = observers.x.au.astype(np.double)
    oy = observers.y.au.astype(np.double)
    oz = observers.z.au.astype(np.double)
    ovx = observers.vx.to(u.au/u.day).value.astype(np.double)
    ovy = observers.vy.to(u.au/u.day).value.astype(np.double)
    ovz = observers.vz.to(u.au/u.day).value.astype(np.double)

    rock = clibspacerocks.py_correct_for_ltt_single_observer(N, x, y, z, vx, vy, vz, ox, oy, oz, ovx, ovy, ovz)
    arr = np.ctypeslib.as_array(rock, (6 * N,)).copy()
    clibspacerocks.free_memory(rock)
  
    x, y, z, vx, vy, vz = arr.reshape(N, 6).T

    x = Distance(x, u.au, allow_negative=True)
    y = Distance(y, u.au, allow_negative=True)
    z = Distance(z, u.au, allow_negative=True)

    return x, y, z, vx * u.au/u.day, vy * u.au/u.day, vz * u.au/u.day


clibspacerocks.py_correct_for_ltt_destnosim.argtypes = [ctypes.c_int,
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'), 
                                              ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

clibspacerocks.py_correct_for_ltt_destnosim.restype = ctypes.POINTER(ctypes.c_double)

def correct_for_ltt_destnosim(x, y, z, vx, vy, vz, ox, oy, oz, ovx, ovy, ovz):

    N = len(x)
    rock = clibspacerocks.py_correct_for_ltt_destnosim(N, x, y, z, vx, vy, vz, ox, oy, oz, ovx, ovy, ovz)
    arr = np.ctypeslib.as_array(rock, (3 * N,)).copy()
    clibspacerocks.free_memory(rock)
  
    x, y, z = arr.reshape(N, 3).T
    return x, y, z

