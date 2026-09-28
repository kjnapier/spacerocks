# SPICE validation against CSPICE

These scripts check the native Rust SPICE reader against NAIF's CSPICE (via `spiceypy`).

1. `gen_kernels.py` uses CSPICE's own writers to create synthetic kernels covering every
   supported segment type (SPK 1, 2, 3, 5, 8, 9, 12, 13, 14, 15, 17, 18, 19, 20, 21;
   binary PCK 2, 3, 20), plus a frame kernel (TK frames, body names), a text PCK with an
   IAU rotation model, and a meta-kernel.
2. `compare.py` evaluates states and frame transforms at hundreds of epochs (including
   segment/record/interval boundaries) with both implementations and reports the worst
   relative difference. With `SPACEROCKS_KERNELS` set it also checks real kernels
   (de440s, jwst_pred, polymele, Earth high-precision BPCs) if they are present.

```bash
pip install spiceypy numpy
cargo build --release --example spice_check
python validation/spice/gen_kernels.py
SPACEROCKS_KERNELS=~/data/spice python validation/spice/compare.py
```

Last run (2026-09-26): all checks pass. Most are bit-identical to CSPICE; the rest
(Lagrange/Hermite equal-step types 8/12, two-body types 5/15/17) agree to <6e-14 relative.
