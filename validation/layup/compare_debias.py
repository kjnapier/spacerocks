"""Star-catalog debiasing: spacerocks' `orbfit.debias` vs layup's `debias`, on synthetic bias
tables (JPL's is not needed): random offsets and proper motions for every catalog and pixel.

1. nside 256 (the real table's resolution), given to layup as its in-memory dictionary and to
   spacerocks as its binary cache. This checks the HEALPix lookup and the correction on 200,000
   random positions plus the poles, the polar-cap boundary (|z| = 2/3), and RA 0/360.
2. nside 64, written as a `bias.dat` text file that both codes parse (layup's
   `generate_bias_dict` and `BiasTable::load`).

    python compare_debias.py
"""
import os, struct, tempfile
import numpy as np

from layup.utilities.debiasing import MPC_CATALOGS, debias as layup_debias, generate_bias_dict
from spacerocks import orbfit

rng = np.random.default_rng(2020)
names = list(MPC_CATALOGS)
codes = list(MPC_CATALOGS.values())


def random_table(nside):
    npix = 12 * nside * nside
    t = np.empty((npix, 4 * len(codes)), dtype=np.float32)
    t[:, 0::4] = rng.uniform(-0.5, 0.5, (npix, len(codes)))    # RA cos Dec offset, arcsec
    t[:, 1::4] = rng.uniform(-0.5, 0.5, (npix, len(codes)))    # Dec offset, arcsec
    t[:, 2::4] = rng.uniform(-8.0, 8.0, (npix, len(codes)))    # proper motions, mas/yr
    t[:, 3::4] = rng.uniform(-8.0, 8.0, (npix, len(codes)))
    return np.round(t, 3).astype(np.float32)


def as_dict(t):
    return {c: {k: t[:, 4 * j + i].astype(np.float64) for i, k in enumerate(("ra", "dec", "pm_ra", "pm_dec"))}
            for j, c in enumerate(codes)}


def points(n):
    ra = rng.uniform(0, 360, n)
    dec = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))
    zb = np.degrees(np.arcsin(2 / 3))
    special_dec = [90.0, -90.0, 89.9999, -89.9999, zb, -zb, zb + 1e-9, zb - 1e-9, 0.0, 41.0, 60.0, -75.0]
    special_ra = [0.0, 360 - 1e-10, 45.0, 90.0, 180.0, 270.0, 359.99999]
    sra, sdec = np.meshgrid(special_ra, special_dec)
    ra = np.concatenate([ra, sra.ravel()])
    dec = np.concatenate([dec, sdec.ravel()])
    jd = rng.uniform(2433282.5, 2462502.5, len(ra))            # 1950-2030
    cats = [(names + codes + ["Gaia2", "UCAC5", None, ""])[i % (2 * len(codes) + 4)] for i in range(len(ra))]
    return ra, dec, jd, cats


def compare(label, bias_dict, nside, table_path, n):
    ra, dec, jd, cats = points(n)
    ref = np.array([layup_debias(r, d, j, c, bias_dict, nside=nside) for r, d, j, c in zip(ra, dec, jd, cats)])
    r2, d2 = orbfit.debias(np.radians(ra), np.radians(dec), jd, cats, table=table_path, timescale="tdb")
    dra = (np.degrees(r2) - ref[:, 0] + 180) % 360 - 180
    ddec = np.degrees(d2) - ref[:, 1]
    shift = np.hypot(((ref[:, 0] - ra + 180) % 360 - 180) * np.cos(np.radians(dec)), ref[:, 1] - dec)
    err = np.hypot(dra * np.cos(np.radians(dec)), ddec) * 3.6e6          # mas
    print(f"{label}: {len(ra)} detections, corrections up to {shift.max() * 3600:.2f}\", "
          f"max difference {err.max():.2e} mas, {np.sum(err > 1e-3)} above 1 uas")


with tempfile.TemporaryDirectory() as tmp:
    # 1. nside 256, binary for spacerocks, dictionary for layup
    t = random_table(256)
    with open(os.path.join(tmp, "bias.bin"), "wb") as f:
        f.write(b"SRBIAS1\0" + struct.pack("<QQ", 256, t.shape[1]))
        f.write(t.astype("<f4").tobytes())
    compare("nside 256", as_dict(t), 256, os.path.join(tmp, "bias.dat"), 200000)

    # 2. nside 64, as JPL's text format, parsed by both
    d64 = os.path.join(tmp, "n64")
    os.makedirs(d64)
    t = random_table(64)
    with open(os.path.join(d64, "bias.dat"), "w") as f:
        for i in range(23):
            f.write(f"! synthetic header line {i + 1}\n")
        np.savetxt(f, t, fmt="%.3f")
    compare("nside 64 (text)", generate_bias_dict(d64), 64, os.path.join(d64, "bias.dat"), 50000)
