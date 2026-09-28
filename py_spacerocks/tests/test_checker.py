"""Smoke tests of spacerocks.checker from Python, on real astrometry.

    SPACEROCKS_KERNELS=/path/to/kernels pytest py_spacerocks/tests

Needs latest_leapseconds.tls, de440(s).bsp, sb441-n16.bsp and an Earth orientation file
(earth_*.bpc); skipped without SPACEROCKS_KERNELS. Uses the Bennu detections of
tests/data/bennu_mpc_2010_2025.csv and the MPC orbit of (3666) Holman in
tests/data/mpc_orb_3666.json. Checks against the full MPCORB are in validation/checker/.
"""
import glob
import json
import os

import numpy as np
import pytest

from spacerocks import checker, orbfit
from spacerocks.spice import SpiceKernel

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "..", "..", "tests", "data")
KDIR = os.environ.get("SPACEROCKS_KERNELS")
pytestmark = pytest.mark.skipif(not KDIR, reason="set SPACEROCKS_KERNELS")
ARCSEC = np.pi / 180 / 3600


@pytest.fixture(scope="module")
def kernel():
    k = SpiceKernel()
    planets = "de440.bsp" if os.path.exists(os.path.join(KDIR, "de440.bsp")) else "de440s.bsp"
    for f in ["latest_leapseconds.tls", planets, "sb441-n16.bsp"]:
        k.load(os.path.join(KDIR, f))
    for f in glob.glob(os.path.join(KDIR, "earth_*.bpc")):
        k.load(f)
    return k


@pytest.fixture(scope="module")
def bennu():
    rows = [line.strip().split(",") for line in open(os.path.join(DATA, "bennu_mpc_2010_2025.csv")).readlines()[1:]]
    col = lambda i: np.array([float(r[i]) for r in rows])
    return dict(epoch=col(0), ra=col(1), dec=col(2), sigma_ra=col(3), sigma_dec=col(4), observer=[r[5] for r in rows])


def run(cat, b, kernel, idx, **kw):
    return checker.check(cat, b["ra"][idx], b["dec"][idx], b["epoch"][idx], list(np.array(b["observer"])[idx]), kernel,
                         sigma_ra=b["sigma_ra"][idx], sigma_dec=b["sigma_dec"][idx], timescale="tdb", **kw)


@pytest.fixture(scope="module")
def catalog(bennu, kernel):
    # Bennu fitted to its detections before 2020, and Holman's MPC orbit.
    early = bennu["epoch"] < 2458849.5
    f = orbfit.fit(bennu["ra"][early], bennu["dec"][early], bennu["epoch"][early], list(np.array(bennu["observer"])[early]), kernel,
                   sigma_ra=bennu["sigma_ra"][early], sigma_dec=bennu["sigma_dec"][early], timescale="tdb", robust=True)
    assert f.flag == 0, f
    cat = checker.Catalog.from_fits({"101955": f}, h=[20.2])
    cat.extend(checker.Catalog.from_mpc_orb(os.path.join(DATA, "mpc_orb_3666.json"), kernel))
    return cat


def test_catalog_basics(catalog, tmp_path):
    assert len(catalog) == 2 and catalog.names == ["101955", "3666"]
    assert catalog.has_covariance.all()
    assert catalog.index("3666") == 1
    path = str(tmp_path / "cat.srcat")
    catalog.save(path)
    back = checker.Catalog.load(path)
    assert back.names == catalog.names
    assert np.array_equal(back.states, catalog.states)
    assert catalog.select(np.array([False, True])).names == ["3666"]


def test_later_detections_are_identified(catalog, bennu, kernel):
    late = np.where(bennu["epoch"] >= 2458849.5)[0]
    r = run(catalog, bennu, kernel, late)
    got = np.zeros(len(late), bool)
    ok = (np.array(r["name"]) == "101955") & r["consistent"]
    got[r["detection"][ok]] = True
    # 20 of 23: the misses are 2024-25 detections reported at 0.2" that lie 1-2" off a
    # gravity-only fit to 2011-2019 (Bennu drifts under the Yarkovsky effect).
    assert got.mean() > 0.8, got.mean()
    assert np.all(np.array(r["name"])[r["consistent"]] == "101955")      # never Holman
    assert np.all(r["from_covariance"]) and np.all(np.isfinite(r["mag"]))
    # one value per pair, ready for pandas
    assert len({len(v) for v in r.values()}) == 1


def test_moved_detections_are_not_consistent(catalog, bennu, kernel):
    idx = np.where(bennu["epoch"] >= 2458849.5)[0][:20]
    b = dict(bennu)
    b["dec"] = bennu["dec"] + 30 * ARCSEC
    r = run(catalog, b, kernel, idx, radius=0.0)
    assert len(r["detection"]) == 0
    r = run(catalog, b, kernel, idx, radius=600.0)
    assert len(r["detection"]) > 0 and not r["consistent"].any()
    assert np.allclose(r["ddec"], 30, atol=5)
