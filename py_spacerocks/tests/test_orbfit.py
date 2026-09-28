"""Smoke tests of the Python orbfit API, end to end on real astrometry.

    SPACEROCKS_KERNELS=/path/to/kernels pytest py_spacerocks/tests

The kernel directory needs latest_leapseconds.tls, de440.bsp (or de440s.bsp), sb441-n16.bsp,
an Earth orientation file (earth_*.bpc) and pck00010.tpc; the tests are skipped without it. The
astrometry is tests/data/bennu_mpc_2010_2025.csv: 310 MPC detections of (101955) Bennu, whose
first opposition alone defeats layup's pipeline.
"""
import glob
import os

import numpy as np
import pytest

from spacerocks import SpaceRock, orbfit
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "..", "..", "tests", "data", "bennu_mpc_2010_2025.csv")
KDIR = os.environ.get("SPACEROCKS_KERNELS")
pytestmark = pytest.mark.skipif(not KDIR, reason="set SPACEROCKS_KERNELS")
ARCSEC = np.pi / 180 / 3600


@pytest.fixture(scope="module")
def kernel():
    k = SpiceKernel()
    planets = "de440.bsp" if os.path.exists(os.path.join(KDIR, "de440.bsp")) else "de440s.bsp"
    for f in ["latest_leapseconds.tls", planets, "sb441-n16.bsp", "pck00010.tpc"]:
        k.load(os.path.join(KDIR, f))
    for f in glob.glob(os.path.join(KDIR, "earth_*.bpc")):
        k.load(f)
    return k


@pytest.fixture(scope="module")
def bennu():
    rows = [line.strip().split(",") for line in open(DATA).readlines()[1:]]
    col = lambda i: np.array([float(r[i]) for r in rows])
    return dict(epoch=col(0), ra=col(1), dec=col(2), sigma_ra=col(3), sigma_dec=col(4), observer=[r[5] for r in rows])


def fit(b, kernel, idx=slice(None), **kw):
    return orbfit.fit(b["ra"][idx], b["dec"][idx], b["epoch"][idx], list(np.array(b["observer"])[idx]), kernel,
                      sigma_ra=b["sigma_ra"][idx], sigma_dec=b["sigma_dec"][idx], timescale="tdb", **kw)


@pytest.fixture(scope="module")
def full(bennu, kernel):
    return fit(bennu, kernel, robust=True)


def test_robust_fit_recovers_bennu(full, bennu):
    assert full.flag == 0, full
    assert full.used.sum() >= 0.95 * len(bennu["ra"])
    rock = full.rock
    assert abs(rock.a() - 1.126) < 0.01                      # barycentric; heliocentric a = 1.1260
    res = full.residuals[full.used, :2]
    assert np.sqrt(np.mean(res ** 2)) < 1.0 * ARCSEC


def test_bk_engine_and_herget_agree_with_the_default(full, bennu, kernel):
    # the same detections (all of them) from the same start, Cartesian vs Bernstein-Khushalani
    cart = fit(bennu, kernel, initial=full.rock)
    bk = fit(bennu, kernel, initial=full.rock, engine="bk_native")
    assert cart.flag == 0 and bk.flag == 0, (cart, bk)
    d = bk.state - cart.state
    assert np.sqrt(d @ np.linalg.solve(cart.state_covariance, d)) < 0.01
    # Herget instead of Gauss for the initial orbit, robust: the same orbit
    h = fit(bennu, kernel, robust=True, iod="herget")
    assert h.flag == 0, h
    rock = h.rock                     # Herget's epoch is its first detection: move to the other's
    rock.propagate(full.epoch, kernel)
    d = np.concatenate([rock.position, rock.velocity]) - full.state
    assert np.sqrt(d @ np.linalg.solve(full.state_covariance, d)) < 0.1
    seed = orbfit.herget_iod(bennu["ra"][:20], bennu["dec"][:20], bennu["epoch"][:20], bennu["observer"][:20], kernel,
                             timescale="tdb")
    assert seed is not None


def test_sequential_update_and_skip(full, bennu, kernel):
    n = len(bennu["ra"])
    old = fit(bennu, kernel, idx=slice(0, n - 40), initial=full.rock)
    assert old.flag == 0
    upd = fit(bennu, kernel, prior=old)
    assert upd.route == "sequential"
    ref = fit(bennu, kernel, initial=old.rock)
    d = upd.state - ref.state
    assert np.sqrt(d @ np.linalg.solve(ref.state_covariance, d)) < 0.05
    same = fit(bennu, kernel, prior=upd)
    assert same.route == "skip" and np.array_equal(same.state, upd.state)
    assert upd.fingerprint == same.fingerprint and len(upd.keys) == n


def test_fit_many_updates_a_catalog(bennu, kernel):
    n = len(bennu["ra"])
    ids = np.where(np.arange(n) < n // 2, "first", "second")
    args = (bennu["ra"], bennu["dec"], bennu["epoch"], bennu["observer"], kernel)
    kw = dict(sigma_ra=bennu["sigma_ra"], sigma_dec=bennu["sigma_dec"], timescale="tdb", robust=True)
    cat = orbfit.fit_many(ids, *args, **kw)
    assert list(cat["id"]) == ["first", "second"]
    again = orbfit.fit_many(ids, *args, prior=cat, **kw)
    assert [r for r, f in zip(again["route"], cat["flag"]) if f == 0] == ["skip"] * int((cat["flag"] == 0).sum())


def test_predict_matches_the_detections(full, bennu, kernel):
    p = orbfit.predict(full, bennu["epoch"], bennu["observer"], kernel, timescale="tdb")
    dra = (p["ra"] - bennu["ra"] + np.pi) % (2 * np.pi) - np.pi
    sep = np.hypot(dra * np.cos(bennu["dec"]), p["dec"] - bennu["dec"])
    assert np.median(sep[full.used]) < 0.5 * ARCSEC
    assert np.all(p["sigma_major"] >= p["sigma_minor"]) and np.all(p["sigma_minor"] > 0)
    later = orbfit.predict(full, [Time(2462000.5, "tdb", "jd")], "X05", kernel)
    assert later["sigma_major"][0] > 0 and later["delta"][0] > 0


def test_comet_orbits(full, kernel):
    assert orbfit.comet_orbits(full, kernel) == {"original": None, "future": None}
    mu = 2.9630927487993194e-4
    r = np.array([1.2, -2.5, 1.1])
    v = np.array([-0.2, 0.9, 0.35]) - r / np.linalg.norm(r) * 0.55
    v *= np.sqrt(2 * mu / np.linalg.norm(r)) * 0.9995 / np.linalg.norm(v)
    rock = SpaceRock.from_xyz("comet", *r, *v, Time(2460400.5, "tdb", "jd"), "J2000", "SSB")
    o = orbfit.comet_orbits(rock, kernel)
    assert o["original"]["reached"] and abs(o["original"]["distance"] - 250) < 1e-3
    assert o["future"]["epoch"] > 2460400.5 > o["original"]["epoch"]


def test_utc_before_1972():
    # 1965 Jan 1: TT - UTC = 32.184 + 3.5401300 s (ERFA); 1900: Delta T = -1.977 s
    for jd, dt in [(2438761.5, 32.184 + 3.54013), (2415020.5, -1.9754), (2451545.0, 64.184)]:
        tt = Time(jd, "utc", "jd").tt().jd()
        assert abs((tt - jd) * 86400 - dt) < 1e-3


def test_weights_observers_occultations(kernel):
    s = orbfit.veres_sigma(["F51", "703", "500"], [Time(2459000.5, "utc", "jd")] * 3)
    assert s.shape == (3,) and np.all(s > 0)
    obs = orbfit.observers(["X05", "C51"], [2460000.5, 2460000.5], kernel, timescale="tdb",
                           sys=[None, "ICRF_KM"], pos=np.array([[np.nan] * 3, [7000.0, 0, 0]]))
    assert obs.shape == (2, 6)
    ra, dec = orbfit.occultation_radec(np.array([1.0]), np.array([0.3]), np.array([0.0]), np.array([0.0]))
    assert abs(ra[0] - 1.0) < 1e-15 and abs(dec[0] - 0.3) < 1e-15
