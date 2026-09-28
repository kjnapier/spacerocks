# Time scales against ERFA

`compare_erfa.py` checks `spacerocks.time` UTC → TT against ERFA (`pip install pyerfa`) from
1960 to 2030, round trips from 1600 to 2030, and TT − UT before 1960 against Stephenson,
Morrison & Hohenkerk's Table S15.2020 (`pip install skyfield`, which bundles it).

Results (2026-09-27):

```
1. TT - UTC, 1960-2030, 6916 epochs: max difference from ERFA 33.8 us (1960-1972: 33.8 us); one ulp of a JD is 40.2 us
2. round trips 1600-2030, 20000 epochs: UTC->TT->UTC max 0.0 us, UTC->TDB->UTC max 0.0 us
3. TT - UT before 1960, 20000 epochs: max difference from Table S15.2020 33.8 us
```

Every difference is within the resolution of a Julian date held in one double. Days that end in
a step of TAI − UTC are left out of (1): ERFA reads a UTC Julian date as a quasi-JD that spreads
the step over that day, while spacerocks applies it at 0h, when it happened.
