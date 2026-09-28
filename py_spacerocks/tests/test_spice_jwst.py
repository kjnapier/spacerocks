"""JWST states from spacerocks compared live against spiceypy (CSPICE).

    SPACEROCKS_JWST_BSP=/path/to/jwst_pred.bsp SPACEROCKS_KERNELS=/path/to/kernels pytest py_spacerocks/tests

Needs spiceypy. The JWST-relative-to-Earth checks need only jwst_pred.bsp (from SPACEROCKS_JWST_BSP,
or jwst_pred.bsp inside SPACEROCKS_KERNELS). The Sun and solar-system-barycentre checks also need
de440s.bsp (or de440.bsp) in SPACEROCKS_KERNELS and are skipped without it.
"""
import os

import numpy as np
import pytest

sp = pytest.importorskip("spiceypy")

from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

KDIR = os.environ.get("SPACEROCKS_KERNELS")
JWST = os.environ.get("SPACEROCKS_JWST_BSP") or (KDIR and os.path.join(KDIR, "jwst_pred.bsp"))
pytestmark = pytest.mark.skipif(not (JWST and os.path.exists(JWST)), reason="set SPACEROCKS_JWST_BSP")


def planets():
    for name in ["de440.bsp", "de440s.bsp"]:
        if KDIR and os.path.exists(os.path.join(KDIR, name)):
            return os.path.join(KDIR, name)
    return None


@pytest.fixture(scope="module")
def cspice():
    sp.kclear()
    sp.furnsh(JWST)
    if planets():
        sp.furnsh(planets())
    yield
    sp.kclear()


@pytest.fixture(scope="module")
def kernel():
    k = SpiceKernel()
    if planets():
        k.load(planets())
    k.load(JWST)  # generic loader must detect the SPK
    return k


def segments():
    h = sp.dafopr(JWST)
    sp.dafbfs(h)
    out = []
    while sp.daffna():
        dc, ic = sp.dafus(sp.dafgs(), 2, 6)
        out.append((dc[0], dc[1], ic[1], ic[3]))
    sp.dafcls(h)
    return out


def epochs(n, lo, hi):
    """Random epochs plus every segment's ends, just inside them, and midpoint."""
    ets = list(np.random.default_rng(0).uniform(lo, hi, n))
    for a, b, *_ in segments():
        if lo <= a and b <= hi:
            ets += [a, b, a + 1e-3, b - 1e-3, 0.5 * (a + b)]
    return ets


def compare(kernel, observer, ets):
    worst_r = worst_v = 0.0
    et_r = et_v = None
    for et in ets:
        t = Time(2451545.0 + et / 86400.0, "tdb", "jd")
        et_used = (t.tdb().jd() - 2451545.0) * 86400.0  # the ET spacerocks actually evaluates
        ours = np.asarray(kernel.state(-170, observer, t, units="km"))
        ref, _ = sp.spkgeo(-170, et_used, "J2000", observer)
        dr, dv = np.linalg.norm(ours[:3] - ref[:3]), np.linalg.norm(ours[3:] - ref[3:])
        if dr > worst_r:
            worst_r, et_r = dr, et
        if dv > worst_v:
            worst_v, et_v = dv, et
    return worst_r, worst_v, et_r, et_v


# Tolerances. On x86-64 spacerocks and CSPICE agree bit for bit. On Apple Silicon CSPICE is
# built with fused multiply-add contraction and Rust is not, and the type-13 Hermite derivative
# amplifies those last-bit differences (observed ~2e-11 km/s). 1 mm and 1 um/s are far below
# anything physical and many orders of magnitude below what a real reader bug produces.
TOL_R = 1e-6   # km
TOL_V = 1e-9   # km/s


def test_file_layout(kernel):
    assert kernel.loaded_kernels[-1][1] == "SPK"  # JWST was loaded last (highest priority)
    segs = segments()
    assert {s[2] for s in segs} == {399} and {s[3] for s in segs} == {13}


def test_relative_to_earth(kernel, cspice):
    segs = segments()
    lo, hi = min(s[0] for s in segs), max(s[1] for s in segs)
    dr, dv, et_r, et_v = compare(kernel, 399, epochs(5000, lo, hi))
    assert dr < TOL_R, f"max |dr| = {dr:.3e} km at ET {et_r!r}"
    assert dv < TOL_V, f"max |dv| = {dv:.3e} km/s at ET {et_v!r}"


@pytest.mark.parametrize("observer", [10, 0], ids=["sun", "ssb"])
def test_relative_to_sun_and_ssb(kernel, cspice, observer):
    if not planets():
        pytest.skip("needs de440s.bsp or de440.bsp in SPACEROCKS_KERNELS")
    segs = segments()
    lo = max(min(s[0] for s in segs), sp.wnfetd(sp.spkcov(planets(), 399), 0)[0])
    hi = min(max(s[1] for s in segs), sp.wnfetd(sp.spkcov(planets(), 399), 0)[1])
    dr, dv, et_r, et_v = compare(kernel, observer, epochs(2000, lo, hi))
    assert dr < 10 * TOL_R, f"max |dr| = {dr:.3e} km at ET {et_r!r}"   # chained through DE
    assert dv < 10 * TOL_V, f"max |dv| = {dv:.3e} km/s at ET {et_v!r}"


def test_outside_coverage_raises(kernel):
    lo = min(s[0] for s in segments())
    with pytest.raises(ValueError):
        kernel.state(-170, 399, Time(2451545.0 + (lo - 86400.0) / 86400.0, "tdb", "jd"))


def test_load_bpc_rejects_spk():
    with pytest.raises(ValueError):
        SpiceKernel().load_bpc(JWST)
