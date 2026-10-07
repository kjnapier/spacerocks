"""RockCollection: shared frame and origin, columns and one-pass elements (no kernels needed)."""
import math

import numpy as np
import pytest

from spacerocks import RockCollection, SpaceRock
from spacerocks.time import Time


def make_rocks(n, epoch, plane="ECLIPJ2000", origin="SUN"):
    rocks = []
    for i in range(n):
        e = 0.0 if i % 5 == 0 else abs(math.sin(0.13 * i)) * 0.95
        rocks.append(SpaceRock.from_kepler(f"r{i}", 1.0 + i, e, 0.3 * i % 3.0, 0.7 * i % 6.2,
                                           1.1 * i % 6.2, 0.9 * i % 6.2 - 3.1, epoch, plane, origin))
    return rocks


@pytest.fixture
def epoch():
    return Time(2460000.5, "tdb", "jd")


def test_columns_and_elements_match_rocks(epoch):
    rocks = make_rocks(30, epoch)
    c = RockCollection()
    for r in rocks:
        c.add(r)
    assert len(c) == 30
    assert c.reference_plane == "ECLIPJ2000"
    assert c.origin == "SUN"

    states = c.states
    assert states.shape == (30, 6)
    assert np.array_equal(states[:, 0], c.x)
    assert np.array_equal(states[:, 5], c.vz)
    assert np.array_equal(c.x, [r.x for r in rocks])

    el = c.elements()
    for key in ["a", "e", "q", "inc", "node", "arg", "true_anomaly", "mean_anomaly"]:
        expected = np.array([getattr(r, key)() for r in rocks])
        assert np.array_equal(el[key], expected), key
        assert np.array_equal(getattr(c, key)(), expected), key
    assert rocks[3].elements()["a"] == rocks[3].a()

    assert c[3].name == "r3"
    assert c[-1].name == "r29"
    assert c.get("r7").x == rocks[7].x
    assert c.name[:2] == ["r0", "r1"]


def test_frames_are_normalised_and_origins_checked(epoch):
    c = RockCollection(reference_plane="J2000", origin="SUN")
    r = make_rocks(1, epoch)[0]
    c.add(r)
    assert c[0].reference_plane == "J2000"
    r2 = SpaceRock.from_xyz("r", r.x, r.y, r.z, r.vx, r.vy, r.vz, epoch, "ECLIPJ2000", "SUN")
    r2.change_reference_plane("J2000")
    assert c[0].x == r2.x

    with pytest.raises(ValueError):
        c.add(make_rocks(1, epoch, origin="SSB")[0])
    assert len(c) == 1


def test_filter_and_propagate(epoch):
    rocks = make_rocks(20, epoch)
    c = RockCollection()
    for r in rocks:
        c.add(r)
    f = c.filter([i % 2 == 0 for i in range(20)])
    assert len(f) == 10

    t1 = Time(2460100.5, "tdb", "jd")
    c.analytic_propagate(t1)
    for i, r in enumerate(rocks):
        r.analytic_propagate(t1)
        assert c[i].x == r.x
    c.change_reference_plane("J2000")
    assert c.reference_plane == "J2000"
