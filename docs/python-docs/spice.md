# SPICE Module

`spacerocks.spice.SpiceKernel` reads NAIF SPICE kernels natively (no CSPICE installation needed).

## What is supported

| Kernel | Support |
|---|---|
| SPK (`.bsp`) | Segment types 1, 2, 3, 5, 8, 9, 12, 13, 14, 15, 17, 18, 19, 20, 21 |
| Binary PCK (`.bpc`) | Types 2, 3, 20 (e.g. `earth_latest_high_prec.bpc`, `earth_*_combined.bpc`) |
| Text kernels | LSK, text PCK (GM values, IAU rotation models), frame kernels (TK frames, body names) |
| Meta-kernels (`.tm`) | `KERNELS_TO_LOAD`, `PATH_SYMBOLS` / `PATH_VALUES`, `+` continuation |
| Frames | 21 built-in inertial frames, `ITRF93` and `IAU_*` body-fixed frames, TK frames |

Both little- and big-endian binary files are read. Type 10 (TLE) SPK segments, CK and DSK
kernels are not supported.

SPICE's precedence rules apply: files loaded later take priority over earlier ones, and within
a file later segments take priority. States are chained between any target and observer using
whichever segments cover the requested epoch, so, for example, JWST relative to the Sun works
with `jwst_pred.bsp` (JWST relative to Earth) plus `de440s.bsp`.

Every `SpiceKernel` is independent. `copy()` makes a new kernel that shares the files already
loaded (memory-mapped, so no extra memory), and you can then load more into either one.

## Default kernels and automatic download

```python
k = SpiceKernel.defaults()                 # download anything missing, then load
k = SpiceKernel.defaults(download=False)   # only use files already on disk
k = SpiceKernel.defaults(update=True)      # also fetch a newer Earth orientation file if available
k = SpiceKernel.from_config("config.toml") # your own kernel set (see config.toml in the repo)
SpiceKernel.cache_dir()                    # where downloads go
```

The default set is `latest_leapseconds.tls`, `de440s.bsp` (planets, 1849–2150), the newest
`earth_1962_*_combined.bpc` (high-precision Earth orientation, history and prediction),
`gm_de440.tpc`, `pck00010.tpc` (the IAU rotation model of the Earth, used for observatory
positions before 1962, where the binary Earth orientation starts), and `sb441-n16.bsp` (the 16
asteroid perturbers used by `propagate`, about 650 MB). Files are cached in `~/.spacerocks/spice`; set `SPACEROCKS_SPICE_DIR` to put the cache
elsewhere (for example on scratch storage on a cluster). Downloads are written to a `.part` file
and renamed when complete, so an interrupted download is never mistaken for a kernel.

## Example

```python
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

k = SpiceKernel()
k.load("~/data/spice/de440s.bsp")          # or a meta-kernel: k.load("my_kernels.tm")
k.load("~/data/spice/earth_latest_high_prec.bpc")
k.load("~/data/spice/jwst_pred.bsp")
print(k)

t = Time(2460500.5, "tdb", "jd")
k.state("JWST", "SUN", t)                       # AU, AU/day in J2000
k.state(-170, 399, t, frame="ECLIPJ2000", units="km")
k.pxform("ITRF93", "J2000", t)                  # 3x3 rotation
k.sxform("ITRF93", "J2000", t)                  # 6x6 state transformation
k.coverage("JWST")                              # [(start_jd, end_jd), ...]

k2 = k.copy()
k2.unload("~/data/spice/jwst_pred.bsp")         # k is unaffected
```

## Methods

| Method | Description |
|---|---|
| `SpiceKernel()` | Empty kernel |
| `SpiceKernel.defaults(download=True, update=False)` | Standard kernel set, downloading missing files |
| `SpiceKernel.from_config(path, download=None)` | Kernel set from a TOML file |
| `load_defaults(download=True, update=False)` | Add the standard set to an existing kernel |
| `SpiceKernel.cache_dir()` | Download/cache directory |
| `load(path)` | Load any supported kernel; the type is detected from the file |
| `load_spk(path)`, `load_pck(path)`, `load_bpc(path)` | Load a specific binary kernel type |
| `unload(path) -> bool` | Unload a file (a meta-kernel also unloads what it loaded) |
| `clear()` | Unload everything |
| `copy()` | New kernel sharing the loaded files |
| `loaded_kernels` | List of `(path, kind)`, lowest priority first |
| `bodies` | NAIF IDs with SPK data |
| `body_id(name)`, `body_name(id)` | Name ↔ ID (built-in NAIF table plus names from loaded text kernels) |
| `coverage(body)` | SPK coverage as TDB Julian date intervals |
| `state(target, observer, epoch, frame="J2000", units="au")` | Geometric state (no light-time correction) |
| `pxform(from, to, epoch)` | Rotation matrix |
| `sxform(from, to, epoch)` | State transformation matrix |

Bodies can be given as names or NAIF IDs. Errors (missing coverage, unknown bodies or frames)
raise `ValueError` with a description of what is missing.
