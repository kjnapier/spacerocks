"""nbody Simulation on its state arrays (no kernels needed)."""
import math

import numpy as np

from spacerocks import SpaceRock
from spacerocks.nbody import Integrator, Simulation
from spacerocks.time import Time

GM_SUN = 0.01720209895 ** 2


def make_sim(integrator):
    t = Time(2460000.5, "tdb", "jd")
    sim = Simulation()
    sim.set_epoch(t)
    sim.set_integrator(integrator)
    sun = SpaceRock.from_xyz("sun", 0, 0, 0, 0, 0, 0, t, "J2000", "SSB")
    sun.set_mass(1.0)
    jup = SpaceRock.from_xyz("jupiter", 5.2, 0, 0, 0, math.sqrt(GM_SUN / 5.2), 0, t, "J2000", "SSB")
    jup.set_mass(9.547919e-4)
    tp = SpaceRock.from_xyz("tp", 0, 30.0, 0, -math.sqrt(GM_SUN / 30.0), 0, 0.001, t, "J2000", "SSB")
    for r in (tp, jup, sun):
        sim.add(r)
    return sim


def test_arrays_and_particles_agree():
    sim = make_sim(Integrator.wisdom_holman(20.0))
    sim.move_to_center_of_mass()
    e0 = sim.energy()
    sim.steps(500)
    assert len(sim) == 3
    assert sim.names == ["sun", "jupiter", "tp"]
    assert np.array_equal(sim.masses, [1.0, 9.547919e-4, 0.0])
    states = sim.states
    assert states.shape == (3, 6)
    p = sim.get_particle("tp")
    assert np.array_equal(states[2], [p.x, p.y, p.z, p.vx, p.vy, p.vz])
    assert p.epoch.jd() == sim.epoch.jd()
    assert abs((sim.energy() - e0) / e0) < 1e-6
    assert [r.name for r in sim.particles()] == sim.names
