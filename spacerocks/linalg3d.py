from numba import jit
import numpy as np
# Optimized pure-Python vector algebra in three dimensions.
# for some reason, jit destroys astropy units. I don't think the speedup is worth the sacrifice.
# jit
def cross(a, b):
    x1 = a[1] * b[2] - a[2] * b[1]
    x2 = a[2] * b[0] - a[0] * b[2]
    x3 = a[0] * b[1] - a[1] * b[0]
    return x1, x2, x3

# @jit
def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

# @jit
def norm(a):
    return (a[0]**2 + a[1]**2 + a[2]**2)**0.5

@jit
def euler_rotation(a, b, c, x):
    # Currently assumes z=0. I'm going to modify this to the genreal case.
    xrot = x[0] * (np.cos(a)*np.cos(c) - np.sin(a)*np.sin(c)*np.cos(b)) - \
           x[1] * (np.sin(a)*np.cos(c) + np.cos(a)*np.sin(c)*np.cos(b))
    yrot = x[0] * (np.cos(a)*np.sin(c) + np.sin(a)*np.cos(c)*np.cos(b)) + \
           x[1] * (np.cos(a)*np.cos(c)*np.cos(b) - np.sin(a)*np.sin(c))
    zrot = x[0] * (np.sin(a)*np.sin(b)) + x[1] * (np.cos(a)*np.sin(b))
    return xrot, yrot, zrot
