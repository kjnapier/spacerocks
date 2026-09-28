"""Vereš et al. (2017) uncertainties: spacerocks' `orbfit.veres_sigma` vs layup's
`astrometric_uncertainty_Veres2017`, over every station the model names (plus others), every
catalog it branches on (as ADES names and MPC codes, plus others and none), several program codes,
and dates on both sides of each date split.

    python compare_veres.py
"""
import itertools
import numpy as np

from layup.utilities.astrometric_uncertainty import astrometric_uncertainty_Veres2017 as layup_sigma
from spacerocks import orbfit

stations = ["703", "691", "644", "704", "G96", "F51", "G45", "699", "D29", "C51", "E12", "608", "J75", "645",
            "673", "689", "950", "H01", "J04", "W84", "G83", "K92", "K93", "Q63", "Q64", "V37", "W85", "W86",
            "W87", "K91", "E10", "F65", "Y28", "568", "T09", "T12", "T14", "309", "X05", "I41", "500", "000"]
catalogs = [None, "", "UCAC4", "PPMXL", "Gaia1", "Gaia2", "Gaia3", "Gaia3E", "USNOB1", "USNOB2", "UCAC2",
            "o", "s", "q", "t", "U", "V", "W", "X", "c", "UNK"]
programs = [None, "", "2", "&", "1", "02"]
splits = (2456658.5, 2452640.5, 2452883.5)
jds = sorted({j + d for j in splits for d in (-1.0, 0.0, 1e-6, 1.0)} | {2440000.5, 2461000.5})

rows = list(itertools.product(stations, catalogs, programs, jds))
st, cat, prog, jd = (list(c) for c in zip(*rows))
mine = orbfit.veres_sigma(st, np.array(jd), catalog=cat, program=prog, timescale="tdb") * 206264.80624709636
theirs = np.array([layup_sigma(s, j, catalog=c, program=p) for s, c, p, j in rows])
bad = np.flatnonzero(np.abs(mine - theirs) > 1e-12)
print(f"{len(rows)} combinations, {len(bad)} differ")
for i in bad[:20]:
    print("  ", rows[i], mine[i], theirs[i])
