"""The orbit-fit dashboard from Python, without kernels (its fits fail, but the page and the
handle work). Fits of real astrometry on the dashboard are in test_orbfit.py."""
import json
import time
import urllib.request

import numpy as np

from spacerocks import orbfit
from spacerocks.spice import SpiceKernel


def post(url, query):
    req = urllib.request.Request(url + "control?" + query, method="POST", headers={"X-Spacerocks": "1"})
    return urllib.request.urlopen(req).read()


def settle(d):
    for _ in range(500):
        if not d.fitting:
            return
        time.sleep(0.01)
    raise AssertionError("still fitting")


def test_dashboard_handle():
    n = 5
    d = orbfit.dashboard(np.linspace(1.0, 1.05, n), np.full(n, 0.1), 2460000.5 + np.arange(n), np.tile([1.0, 0.0, 0.0], (n, 1)),
                         SpiceKernel(), timescale="tdb", asteroids=False, name="test rock", open_browser=False)
    assert d.url.startswith("http://127.0.0.1:")
    settle(d)
    assert d.fit is not None and not d.fit.converged
    state = json.loads(urllib.request.urlopen(d.url + "state.json").read())
    assert state["name"] == "test rock" and state["settings"]["asteroids"] is False

    post(d.url, "action=exclude&index=0&on=1")
    assert d.excluded == [True, False, False, False, False]
    post(d.url, "action=done")
    settle(d)
    assert d.wait() is not None
    d.close()
    assert "closed" in repr(d)


def test_dashboard_load_and_command(capfd):
    n = 5
    d = orbfit.dashboard(np.linspace(1.0, 1.05, n), np.full(n, 0.1), 2460000.5 + np.arange(n), np.tile([1.0, 0.0, 0.0], (n, 1)),
                         SpiceKernel(), timescale="tdb", asteroids=False, name="test rock", open_browser=False)
    settle(d)
    assert d.name == "test rock" and not d.loading
    d.close()

    orbfit.main(["--help"])
    assert "spacerocks-fit" in capfd.readouterr().out
    try:
        orbfit.main(["--frobnicate"])
    except SystemExit as e:
        assert e.code == 1
    else:
        raise AssertionError("an unknown option should exit")
