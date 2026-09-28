# SPICE (Rust)

`spacerocks::spice` is a native Rust implementation of the parts of the SPICE toolkit that
spacerocks needs. See the module documentation (`cargo doc --open`) for details.

## Layout

| File | Contents |
|---|---|
| `daf.rs` | Generic DAF reader: memory-mapped, zero-copy for native byte order, byte-swapped otherwise |
| `spk.rs` | SPK segment parsing and evaluation (types 1, 2, 3, 5, 8, 9, 12, 13, 14, 15, 17, 18, 19, 20, 21) |
| `pck.rs` | Binary PCK segments (types 2, 3, 20) |
| `generic.rs` | NAIF generic-segment layout (SPK 14, PCK 3) |
| `text.rs` | Text-kernel parser and `KernelPool` |
| `frames.rs` | Built-in inertial frames, PCK frames (binary and IAU text models), TK frames, `sxform`/`pxform` |
| `bodies.rs` | NAIF body names/IDs (built-in table plus kernel-defined names), `SpiceBody` |
| `config.rs` | Kernel sets, TOML configuration, cache directory and downloading |
| `kernel.rs` | `SpiceKernel`: loading, precedence, chaining, batch queries |
| `math.rs` | Chebyshev, Hermite/Lagrange interpolation, two-body propagation (ports of the SPICELIB routines) |

## Key API

```rust
use spacerocks::spice::{SpiceKernel, et_from_jd};

let mut k = SpiceKernel::new();
k.load("de440s.bsp")?;              // SPK, PCK, text or meta-kernel
k.load("earth_latest_high_prec.bpc")?;

let et = et_from_jd(2_460_500.5);   // TDB seconds past J2000
let s = k.spkgeo(399, 0, et)?;      // km, km/s, J2000 (SPICE spkgeo)
let s = k.state(-170, 10, "ECLIPJ2000", et)?;
let r = k.pxform("ITRF93", "J2000", et)?;   // nalgebra::Matrix3
let x = k.sxform("ITRF93", "J2000", et)?;   // StateTransform { rotation, rate }

// spacerocks units (AU, AU/day, TDB Julian date)
let s = k.state_au(399, 0, 2_460_500.5)?;
let mut out = vec![[0.0; 6]; 3];
k.barycentric_states_au(&[10, 399, 301], 2_460_500.5, &mut out)?;
```

## Default kernel set and downloads

```rust
let k = SpiceKernel::defaults()?;                      // downloads missing files into the cache
let k = SpiceKernel::defaults_with_download(false)?;   // local files only
let k = SpiceKernel::from_config_file("config.toml")?; // custom set; see `spice::config`
```

`spice::config::Config` holds the kernel list (`KernelSpec`s: NAIF subdirectory, explicit URL,
or local-only; `*` selects the newest match), search paths, the cache directory
(`$SPACEROCKS_SPICE_DIR` or `~/.spacerocks/spice`), `auto_download`, `check_for_updates`, and
a NAIF mirror URL.

## Design notes

* **Independent contexts.** There is no global kernel pool. A `SpiceKernel` is immutable while
  queried, so it is `Send + Sync` and can be shared across threads with no locking. Loaded files
  are `Arc`-shared, so `clone()` is cheap and clones can diverge.
* **SPICE semantics.** Later files and later segments take precedence. Coverage is checked per
  segment at the requested epoch, and target/observer chains are resolved the way `SPKGEO`
  resolves them.
* **Numerics.** Evaluation routines follow the SPICELIB algorithms, including their operation
  order where it matters (`CHBINT`, `HRMINT`, `PROP2B`, `TISBOD`, f2c `MOD`). Results are
  bit-identical to CSPICE for most segment types and within 1e-13 relative for the rest; see
  `validation/spice`.
* **Robustness.** Malformed files produce `SpiceError`s, not panics. This was fuzz-tested with
  100k mutated kernels.

## Not supported (yet)

SPK type 10 (TLEs) and 16, CK and DSK kernels, dynamic frames, and aberration corrections in the
SPICE module itself (light-time correction lives in `transforms::correct_for_ltt`).
