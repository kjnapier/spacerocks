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


def test_monitor_page_reports_progress():
    import json
    import urllib.request

    sim = make_sim(Integrator.ias15(1.0))
    url = sim.monitor(open_browser=False)
    assert sim.monitor_url == url
    target = Time(2460000.5 + 1000.0, "tdb", "jd")
    sim.integrate(target)
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    with opener.open(url + "progress.json") as r:
        status = json.load(r)
    assert status["state"] == "done"
    assert status["t_target"] == 1000.0
    assert status["n_particles"] == 3
    assert status["integrator"] == "IAS15"
    sim.stop_monitor()
    assert sim.monitor_url is None


def test_monitor_tracks_elements():
    import json
    import urllib.request

    sim = make_sim(Integrator.ias15(1.0))
    url = sim.monitor(open_browser=False, track="jupiter")
    sim.integrate(Time(2460000.5 + 100.0, "tdb", "jd"))
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    with opener.open(url + "progress.json") as r:
        status = json.load(r)
    assert status["tracked"] == "jupiter"
    assert abs(status["elements"]["a"] - 5.2) < 0.01
    sim.monitor_track(None)
    with opener.open(url + "progress.json") as r:
        assert json.load(r)["tracked"] is None


def _interrupt_after(seconds):
    import os
    import signal
    import threading

    timer = threading.Timer(seconds, lambda: os.kill(os.getpid(), signal.SIGINT))
    timer.start()
    return timer


def test_ctrl_c_stops_integrate():
    import time

    import pytest

    sim = make_sim(Integrator.ias15(1.0))
    target = Time(2460000.5 + 1e9, "tdb", "jd")  # far longer than the test may run
    _interrupt_after(0.3)
    start = time.monotonic()
    with pytest.raises(KeyboardInterrupt):
        sim.integrate(target)
    assert time.monotonic() - start < 5
    assert 2460000.5 < sim.epoch.jd() < target.jd()
    # The simulation is usable afterwards.
    jd = sim.epoch.jd()
    sim.integrate(Time(jd + 10.0, "tdb", "jd"))
    assert abs(sim.epoch.jd() - (jd + 10.0)) < 1e-9


def test_ctrl_c_stops_a_paused_monitor_run():
    import time

    import pytest

    sim = make_sim(Integrator.wisdom_holman(20.0))
    sim.monitor(open_browser=False, paused=True)
    _interrupt_after(0.3)
    start = time.monotonic()
    with pytest.raises(KeyboardInterrupt):
        sim.integrate(Time(2460000.5 + 1e6, "tdb", "jd"))
    assert time.monotonic() - start < 5
    assert sim.epoch.jd() == 2460000.5
    sim.stop_monitor()


def test_run_goes_until_stopped_from_the_page():
    import json
    import threading
    import urllib.request

    sim = make_sim(Integrator.wisdom_holman(1.0))
    url = sim.monitor(open_browser=False)
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))

    def stop():
        req = urllib.request.Request(url + "control?action=stop", method="POST", headers={"X-Spacerocks": "1"})
        opener.open(req).read()

    timer = threading.Timer(0.3, stop)
    timer.start()
    sim.run()
    timer.join()
    assert sim.epoch.jd() > 2460000.5
    with opener.open(url + "progress.json") as r:
        assert json.load(r)["state"] == "stopped"
    sim.stop_monitor()


def test_gr_and_perturber_toggles():
    sim = make_sim(Integrator.ias15(1.0))
    assert sim.solar_gr is False
    sim.set_solar_gr(True)
    assert sim.solar_gr is True
    assert sim.perturbers == [("jupiter", True)]
    sim.set_perturber("jupiter", False)
    assert sim.perturbers == [("jupiter", False)]
    assert sim.masses[1] == 0.0
    assert sim.get_particle("jupiter").mass == 9.547919e-4
    sim.set_perturber("jupiter", True)
    assert sim.masses[1] == 9.547919e-4
