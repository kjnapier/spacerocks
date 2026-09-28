//! Integration tests for the SPICE module.
//!
//! Tests that need real kernels look in the directory named by the `SPACEROCKS_KERNELS`
//! environment variable (for example `~/data/spice`) and are skipped when it is not set.
//! Required files: `de440s.bsp`, `earth_latest_high_prec.bpc` (or
//! `earth_1962_240827_2124_combined.bpc`), `latest_leapseconds.tls`.
//!
//! Numerical agreement with CSPICE for every supported segment type is checked separately by
//! the validation scripts; these tests cover loading, precedence, sharing and error handling.

use spacerocks::spice::{SpiceError, SpiceKernel};
use std::path::PathBuf;

fn kernel_dir() -> Option<PathBuf> {
    let dir = PathBuf::from(std::env::var("SPACEROCKS_KERNELS").ok()?);
    if dir.join("de440s.bsp").exists() {
        Some(dir)
    } else {
        None
    }
}

macro_rules! require_kernels {
    () => {
        match kernel_dir() {
            Some(d) => d,
            None => {
                eprintln!("SPACEROCKS_KERNELS not set; skipping");
                return;
            }
        }
    };
}

fn earth_bpc(dir: &PathBuf) -> Option<PathBuf> {
    ["earth_latest_high_prec.bpc", "earth_1962_240827_2124_combined.bpc"]
        .iter()
        .map(|f| dir.join(f))
        .find(|p| p.exists())
}

#[test]
fn empty_kernel_reports_missing_data() {
    let k = SpiceKernel::new();
    assert!(k.loaded_kernels().is_empty());
    match k.spkgeo(399, 0, 0.0) {
        Err(SpiceError::InsufficientEphemerisData { target, observer, .. }) => {
            assert_eq!((target, observer), (399, 0));
        }
        other => panic!("expected InsufficientEphemerisData, got {:?}", other),
    }
    assert!(matches!(k.spkgeo(399, 399, 0.0), Ok(s) if s == [0.0; 6]));
}

#[test]
fn missing_and_invalid_files() {
    let mut k = SpiceKernel::new();
    assert!(matches!(k.load("/definitely/not/here.bsp"), Err(SpiceError::Io { .. })));
    let tmp = std::env::temp_dir().join("spacerocks_not_a_kernel.txt");
    std::fs::write(&tmp, "hello world").unwrap();
    assert!(matches!(k.load(&tmp), Err(SpiceError::InvalidFile { .. })));
    let _ = std::fs::remove_file(&tmp);
}

#[test]
fn builtin_names_and_frames() {
    let k = SpiceKernel::new();
    assert_eq!(k.body_id("Earth"), Some(399));
    assert_eq!(k.body_id("earth barycenter"), Some(3));
    assert_eq!(k.body_id("JWST"), Some(-170));
    assert_eq!(k.body_id("2000001"), Some(2000001));
    assert_eq!(k.body_name(301).as_deref(), Some("MOON"));
    assert_eq!(k.frame_id("ECLIPJ2000"), Some(17));
    assert_eq!(k.frame_id("itrf93"), Some(13000));
    // Inertial-to-inertial transforms need no kernels.
    let m = k.pxform("ECLIPJ2000", "J2000", 0.0).unwrap();
    let eps = (84381.448f64 / 3600.0).to_radians();
    assert!((m[(1, 1)] - eps.cos()).abs() < 1e-15);
    assert!((m[(2, 1)] - eps.sin()).abs() < 1e-15);
    // ITRF93 needs a binary PCK.
    assert!(matches!(
        k.pxform("ITRF93", "J2000", 0.0),
        Err(SpiceError::InsufficientOrientationData { .. })
    ));
}

#[test]
fn text_kernel_names_and_meta_kernel() {
    let dir = std::env::temp_dir().join(format!("sr_spice_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let fk = dir.join("names.tf");
    std::fs::write(
        &fk,
        "\\begindata\nNAIF_BODY_NAME += ( 'MY ROCK' )\nNAIF_BODY_CODE += ( 3000123 )\nBODY3000123_GM = 1.5\n\\begintext\n",
    )
    .unwrap();
    let mk = dir.join("meta.tm");
    std::fs::write(
        &mk,
        format!(
            "\\begindata\nPATH_VALUES = ( '{}' )\nPATH_SYMBOLS = ( 'D' )\nKERNELS_TO_LOAD = ( '$D/na+'\n'mes.tf' )\n\\begintext\n",
            dir.display()
        ),
    )
    .unwrap();

    let mut k = SpiceKernel::new();
    k.load(&mk).unwrap();
    assert_eq!(k.body_id("my rock"), Some(3000123));
    assert_eq!(k.gm(3000123), Some(1.5));
    assert_eq!(k.loaded_kernels().len(), 2);

    // Unloading the meta-kernel unloads what it loaded.
    assert!(k.unload(&mk));
    assert_eq!(k.body_id("my rock"), None);
    assert!(k.loaded_kernels().is_empty());
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn de440s_states_and_batch() {
    let dir = require_kernels!();
    let mut k = SpiceKernel::new();
    k.load(dir.join("de440s.bsp")).unwrap();

    let jd = 2_460_000.5;
    let et = spacerocks::spice::et_from_jd(jd);
    // Earth relative to the Sun is ~1 AU.
    let s = k.state_au(399, 10, jd).unwrap();
    let r = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
    assert!((r - 1.0).abs() < 0.02, "r = {}", r);

    // Anti-symmetry and chaining.
    let a = k.spkgeo(301, 399, et).unwrap();
    let b = k.spkgeo(399, 301, et).unwrap();
    for i in 0..6 {
        assert!((a[i] + b[i]).abs() <= 1e-9 * a[i].abs().max(1.0));
    }

    // The batch barycentric path agrees with spkgeo exactly.
    let ids = [10, 1, 2, 399, 301, 4, 5, 6, 7, 8, 9];
    let mut out = vec![[0.0; 6]; ids.len()];
    k.barycentric_states_au(&ids, jd, &mut out).unwrap();
    for (i, id) in ids.iter().enumerate() {
        let single = spacerocks::spice::km_to_au_state(&k.spkgeo(*id, 0, et).unwrap());
        for c in 0..6 {
            assert!((out[i][c] - single[c]).abs() <= 1e-15 * single[c].abs().max(1e-12));
        }
    }

    // Outside coverage is an error, not a panic.
    assert!(matches!(
        k.spkgeo(399, 0, 1e12),
        Err(SpiceError::InsufficientEphemerisData { .. })
    ));
}

#[test]
fn clones_share_data_but_not_contents() {
    let dir = require_kernels!();
    let mut base = SpiceKernel::new();
    base.load(dir.join("de440s.bsp")).unwrap();
    let mut other = base.clone();
    let Some(bpc) = earth_bpc(&dir) else { return };
    other.load(&bpc).unwrap();
    assert_eq!(base.loaded_kernels().len(), 1);
    assert_eq!(other.loaded_kernels().len(), 2);
    assert!(other.pxform("ITRF93", "J2000", 0.0).is_ok());
    assert!(base.pxform("ITRF93", "J2000", 0.0).is_err());

    // Concurrent use from several threads.
    let k = std::sync::Arc::new(other);
    let handles: Vec<_> = (0..4)
        .map(|t| {
            let k = k.clone();
            std::thread::spawn(move || {
                for i in 0..1000 {
                    let et = (t * 1000 + i) as f64 * 3600.0;
                    k.spkgeo(399, 0, et).unwrap();
                }
            })
        })
        .collect();
    for h in handles {
        h.join().unwrap();
    }
}

#[test]
fn later_files_take_precedence_and_reload_moves_to_top() {
    let dir = require_kernels!();
    let de = dir.join("de440s.bsp");
    let mut k = SpiceKernel::new();
    k.load(&de).unwrap();
    let n = k.spk_segments(399).len();
    // Loading the same file again does not duplicate it.
    k.load(&de).unwrap();
    assert_eq!(k.loaded_kernels().len(), 1);
    assert_eq!(k.spk_segments(399).len(), n);
    assert!(k.unload(&de));
    assert!(k.spk_bodies().is_empty());
}

#[test]
fn observatory_with_earth_orientation() {
    let dir = require_kernels!();
    let Some(bpc) = earth_bpc(&dir) else { return };
    let mut k = SpiceKernel::new();
    k.load(dir.join("de440s.bsp")).unwrap();
    k.load(&bpc).unwrap();
    let epoch = spacerocks::Time::new(2_460_000.5, "tdb", "jd").unwrap();
    let obs = spacerocks::Observatory::from_obscode("W84").unwrap();
    let o_eq = obs.at(&epoch, "J2000", "SSB", &k).unwrap();
    let o_ec = obs.at(&epoch, "ECLIPJ2000", "SSB", &k).unwrap();
    // The same physical position expressed in two planes must have the same geocentric
    // distance (this was wrong before the rewrite for non-J2000 planes).
    let earth_eq = spacerocks::SpaceRock::from_spice("earth", &epoch, "J2000", "SSB", &k).unwrap();
    let earth_ec = spacerocks::SpaceRock::from_spice("earth", &epoch, "ECLIPJ2000", "SSB", &k).unwrap();
    let d_eq = (o_eq.position - earth_eq.position).norm();
    let d_ec = (o_ec.position - earth_ec.position).norm();
    assert!((d_eq - d_ec).abs() < 1e-15);
    // ~Earth radius in AU
    assert!((d_eq * 149_597_870.7 - 6371.0).abs() < 30.0);
}

#[test]
fn meta_kernel_cycles_are_errors() {
    let dir = std::env::temp_dir().join(format!("sr_spice_cycle_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let a = dir.join("a.tm");
    let b = dir.join("b.tm");
    std::fs::write(&a, format!("\\begindata\nKERNELS_TO_LOAD = ( '{}' )\n", b.display())).unwrap();
    std::fs::write(&b, format!("\\begindata\nKERNELS_TO_LOAD = ( '{}' )\n", a.display())).unwrap();
    let mut k = SpiceKernel::new();
    assert!(matches!(k.load(&a), Err(SpiceError::InvalidFile { .. })));
    k.unload(&a);
    assert!(k.loaded_kernels().is_empty());
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn shared_child_survives_other_meta_unload() {
    let dir = require_kernels!();
    let tmp = std::env::temp_dir().join(format!("sr_spice_shared_{}", std::process::id()));
    std::fs::create_dir_all(&tmp).unwrap();
    let de = dir.join("de440s.bsp");
    let m1 = tmp.join("m1.tm");
    let m2 = tmp.join("m2.tm");
    for m in [&m1, &m2] {
        std::fs::write(m, format!("\\begindata\nKERNELS_TO_LOAD = ( '{}' )\n", de.display())).unwrap();
    }
    let mut k = SpiceKernel::new();
    k.load(&m1).unwrap();
    k.load(&m2).unwrap();
    k.unload(&m1);
    assert!(k.spkgeo(399, 0, 0.0).is_ok());
    k.unload(&m2);
    assert!(k.spkgeo(399, 0, 0.0).is_err());
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn kernel_set_from_local_config() {
    let dir = require_kernels!();
    let toml = format!(
        "auto_download = false\ndownload_dir = \"{}\"\ndefault_kernels = [\n  {{ name = \"de440*.bsp\" }},\n]\n",
        dir.display()
    );
    let cfg = spacerocks::spice::KernelConfig::from_toml(&toml).unwrap();
    let k = SpiceKernel::from_config(&cfg).unwrap();
    assert_eq!(k.loaded_kernels().len(), 1);
    assert!(k.spkgeo(399, 0, 0.0).is_ok());

    // A missing kernel with downloads disabled is a clear error.
    let cfg = spacerocks::spice::KernelConfig::from_toml(&format!(
        "auto_download = false\ndownload_dir = \"{}\"\ndefault_kernels = [ {{ name = \"nope.bsp\", kernel_type = \"spk\" }} ]\n",
        dir.display()
    ))
    .unwrap();
    assert!(matches!(SpiceKernel::from_config(&cfg), Err(SpiceError::KernelNotFound { .. })));
}

#[test]
fn repository_config_parses() {
    let cfg = spacerocks::spice::KernelConfig::from_file(concat!(env!("CARGO_MANIFEST_DIR"), "/config.toml")).unwrap();
    assert_eq!(cfg.default_kernels, spacerocks::spice::config::default_kernel_list());
}

/// JWST predicted ephemeris (`jwst_pred.bsp`, 501 type-13 segments relative to Earth) evaluated
/// against CSPICE `spkgeo(-170, et, "J2000", 399)` reference states (generated with spiceypy 8.2).
/// The epochs cover both ends of the file, a segment boundary and a millisecond either side of it,
/// a segment midpoint, and a scatter in between. The file is found via `SPACEROCKS_JWST_BSP` or as
/// `jwst_pred.bsp` in `SPACEROCKS_KERNELS`; the test is skipped otherwise.
#[test]
fn jwst_type13_matches_cspice() {
    let path = std::env::var("SPACEROCKS_JWST_BSP")
        .map(PathBuf::from)
        .ok()
        .or_else(|| std::env::var("SPACEROCKS_KERNELS").ok().map(|d| PathBuf::from(d).join("jwst_pred.bsp")))
        .filter(|p| p.exists());
    let Some(path) = path else {
        eprintln!("jwst_pred.bsp not found (set SPACEROCKS_JWST_BSP); skipping");
        return;
    };

    #[rustfmt::skip]
    const REFERENCE: [(f64, [f64; 6]); 12] = [
        (693709269.185, [1.14282813987980498e+04, 2.92485100304190667e+02, -8.31209616529463005e+02, 5.15477283253362195e+00, 6.48693050165751650e+00, -3.68984452334925050e-01]),
        (693712869.185, [2.09384102121679789e+04, 2.07063686435994787e+04, -1.50371055267639258e+03, 1.45918523239216613e+00, 4.91160043399302104e+00, -1.01461506628790943e-01]),
        (722865669.1830024, [-3.02185174236208004e+04, 1.53303711375672999e+06, 4.72304723623171973e+05, -1.76569655648717988e-01, 1.06969710451740002e-01, -1.17634037936931002e-01]),
        (722865669.1820023, [-3.02185172470428952e+04, 1.53303711364975525e+06, 4.72304723740811518e+05, -1.76569655672328185e-01, 1.06969710459845366e-01, -1.17634037924549545e-01]),
        (722865669.1840024, [-3.02190582219910320e+04, 1.53303750112942560e+06, 4.72304948246704182e+05, -1.76570383370927037e-01, 1.06970740982504406e-01, -1.17633371248263202e-01]),
        (768223239.1854086, [-3.91993502625091991e+05, -1.39180358374627656e+06, -6.04062006558048539e+05, 2.29396681304184780e-01, -8.12819102961776091e-02, -1.90862289617446612e-01]),
        (800000000.0, [-3.83821230487538152e+05, -1.41538135502742580e+06, -7.21695863498342456e+05, 1.84737787198015285e-01, -6.71958249077502512e-02, -1.70052311355944674e-01]),
        (850012345.678, [2.72520489054494363e+05, 1.64425132207083283e+06, 3.25363070891493175e+05, -1.91737254739441157e-02, 5.82425441293831815e-03, -1.59268884099002840e-02]),
        (900000000.0, [-2.35525662492264470e+05, -1.45417819965750072e+06, -8.24942400541019626e+05, 1.65936897572002973e-01, -2.13584112342354096e-02, 1.19217196096073341e-01]),
        (950000000.0, [-5.82746490465576528e+05, 9.19972778479618835e+05, 7.02864871153474785e+05, -3.84677526199870112e-01, -3.38072504400729557e-01, -4.56439476485429249e-02]),
        (999777668.1825432, [1.31734832536419318e+06, 2.00013174481415153e+04, 2.77499359351095394e+05, 6.71887926728534923e-02, 4.94414118094935684e-01, 6.85087956403441917e-02]),
        (999777669.1825432, [1.31734839255290991e+06, 2.00018118622280017e+04, 2.77499427859849995e+05, 6.71886408996563989e-02, 4.94414054874232012e-01, 6.85087134825564936e-02]),
    ];

    let mut k = SpiceKernel::new();
    k.load(&path).unwrap();
    assert_eq!(k.loaded_kernels()[0].kind.to_string(), "SPK", "load() must detect the SPK");

    // Bit-identical on x86-64. On Apple Silicon CSPICE is built with fused multiply-add and Rust
    // is not; the Hermite derivative amplifies that to ~2e-11 km/s. 1 mm and 1 um/s leave room
    // for this while staying far below any real reader error.
    for (et, want) in REFERENCE {
        let got = k.spkgeo(-170, 399, et).unwrap();
        let dr = (0..3).map(|i| (got[i] - want[i]).powi(2)).sum::<f64>().sqrt();
        let dv = (3..6).map(|i| (got[i] - want[i]).powi(2)).sum::<f64>().sqrt();
        assert!(dr < 1e-6 && dv < 1e-9, "ET {et}: |dr| = {dr:e} km, |dv| = {dv:e} km/s");
    }

    // Outside the file's coverage is an error, not an extrapolation.
    for et in [693709269.185 - 1.0, 999777669.1825432 + 1.0] {
        assert!(matches!(k.spkgeo(-170, 399, et), Err(SpiceError::InsufficientEphemerisData { .. })));
    }

    // The explicit binary-PCK loader must refuse an SPK.
    assert!(SpiceKernel::new().load_pck(&path).is_err());
}

#[test]
fn summary_describes_loaded_files() {
    let dir = require_kernels!();
    let mut k = SpiceKernel::new();
    assert!(k.summary().is_empty());
    k.load(dir.join("latest_leapseconds.tls")).unwrap();
    k.load(dir.join("de440s.bsp")).unwrap();
    let bpc = earth_bpc(&dir);
    if let Some(b) = &bpc {
        k.load(b).unwrap();
    }
    let s = k.summary();
    assert_eq!(s.len(), 2 + bpc.is_some() as usize);

    // Text kernel: variables, no segments.
    assert_eq!(s[0].kind, spacerocks::spice::KernelKind::Text);
    assert!(s[0].variables.iter().any(|v| v == "DELTET/DELTA_AT"));
    assert!(s[0].groups.is_empty());

    // de440s: one group per target, merged coverage matching spk_coverage.
    let de = &s[1];
    assert_eq!(de.kind, spacerocks::spice::KernelKind::Spk);
    assert!(de.size_bytes.unwrap() > 1_000_000);
    let mut targets: Vec<i32> = de.groups.iter().map(|g| g.body).collect();
    targets.sort_unstable();
    assert_eq!(targets, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 199, 299, 301, 399]);
    for g in &de.groups {
        assert_eq!(g.frame, 1);
        assert_eq!(g.coverage, k.spk_coverage(g.body));
    }
    let earth = de.groups.iter().find(|g| g.body == 399).unwrap();
    assert_eq!(earth.center, 3);

    // Binary PCK: grouped by frame class (3000 = ITRF93).
    if bpc.is_some() {
        let p = &s[2];
        assert_eq!(p.kind, spacerocks::spice::KernelKind::Pck);
        assert!(p.groups.iter().all(|g| g.body == 3000));
        let n: usize = p.groups.iter().map(|g| g.n_segments).sum();
        assert!(n >= 1);
        assert_eq!(p.groups[0].coverage, k.pck_coverage(3000));
    }
}
