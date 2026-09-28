"""Run layup's orbit fit on a CSV of astrometry and save the results (reference for spacerocks).

    LAYUP_CACHE=/path/to/layup/cache python run_layup.py input.csv output.npz [n_objects]

Set CONV_FRAC to run with layup's scaled convergence test (default 0, the absolute test), and
WEIGHT_DATA=veres for its Vereš et al. (2017) weights (weight_data=True), and IOD=herget or
IOD=gauss for its other initial-orbit methods (default auto), and ENGINE=bk_native for its
Bernstein-Khushalani fitting engine. With IOD=herget, layup's
herget_iod is wrapped so that a KeplerConvergenceError from its two-body solver fails only that
starting range (unwrapped, it ends the run), restoring the detection times (the error leaves them
light-time shifted), and with the ephemeris from LAYUP_CACHE.
"""
import os, sys, time
import numpy as np
from layup.utilities.file_io.CSVReader import CSVDataReader
from layup.orbitfit import orbitfit
from layup.utilities.data_processing_utilities import LayupObservatory
import spiceypy as spice
import numpy.lib.recfunctions as rfn

def patch_herget(cache):
    import assist
    from layup import iod as layup_iod
    from layup.utilities.herget_iod import herget_with_assist
    from layup.utilities.universal_kepler import KeplerConvergenceError
    ephem = assist.Ephem(os.path.join(cache, "linux_p1550p2650.440"), os.path.join(cache, "sb441-n16.bsp"))

    def herget_iod(observations, seq):
        epochs = [o.epoch for o in observations]
        for rho in (2, 5, 40):
            try:
                s = herget_with_assist(observations, seq, ephem, initial_rho=rho)
            except KeplerConvergenceError:
                for o, e in zip(observations, epochs):
                    o.epoch = e
                continue
            if s:
                return s
        return []
    layup_iod.register_iod("herget", herget_iod)


# At import, so that layup's worker processes (started with "spawn", which re-imports this
# script) get the patch too.
if os.environ.get("IOD") == "herget":
    patch_herget(os.environ["LAYUP_CACHE"])


def main():
    cache = os.environ["LAYUP_CACHE"]
    iod = os.environ.get("IOD", "auto")
    inp, out = sys.argv[1], sys.argv[2]
    nobj = int(sys.argv[3]) if len(sys.argv) > 3 else None

    reader = CSVDataReader(inp, "csv", primary_id_column_name="provID")
    data = reader.read_rows()
    ids = list(dict.fromkeys(data["provID"]))
    if nobj:
        ids = ids[:nobj]
        data = data[np.isin(data["provID"], ids)]

    # observer states exactly as layup computes them (also saved for the comparison)
    obs = LayupObservatory(cache_dir=cache)
    et = np.array([spice.str2et(t) for t in data["obsTime"]])
    d2 = rfn.append_fields(data, "et", et, usemask=False, asrecarray=True)
    pv = obs.obscodes_to_barycentric(d2)

    rows = []
    t0 = time.time()
    for oid in ids:
        t = time.time()
        res = orbitfit(data[data["provID"] == oid], cache, conv_frac=float(os.environ.get("CONV_FRAC", "0")),
                       weight_data=os.environ.get("WEIGHT_DATA", "") == "veres", iod=iod,
                       engine=os.environ.get("ENGINE", "cartesian"))
        rows.append(res)
        print(oid, int(res["flag"][0]), float(res["csq"][0]), int(res["ndof"][0]), int(res["niter"][0]), f"{time.time()-t:.2f}s", flush=True)
    res = np.concatenate(rows)
    np.savez(out, fits=res, obs=data, et=et, pv=pv)
    print("total", time.time() - t0)


if __name__ == "__main__":
    main()
