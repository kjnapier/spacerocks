//! Fit every object in a detections file with `spacerocks::orbfit` and write the results in
//! the layout of `validation/layup/export_inputs.py`'s `layup_fits.csv`.
//!
//!     SPACEROCKS_KERNELS=/path/to/kernels cargo run --release --example orbfit_layup -- \
//!         detections.csv out.csv [observers]
//!
//! `detections.csv` has columns id, epoch_tdb, ra, dec (radians), ox, oy, oz (barycentric
//! observer, AU), stn. With the optional third argument `spacerocks`, observer positions are
//! recomputed from the station codes with `spacerocks::Observatory` instead of taken from the
//! file. The kernel directory must hold de440.bsp, sb441-n16.bsp and Earth orientation files.

use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use rayon::prelude::*;

use spacerocks::orbfit::{determine_orbit, observer_positions, Astrometry, FitOptions};
use spacerocks::SpiceKernel;

fn kernel() -> SpiceKernel {
    let dir = PathBuf::from(std::env::var("SPACEROCKS_KERNELS").expect("set SPACEROCKS_KERNELS"));
    let mut k = SpiceKernel::new();
    for f in ["de440.bsp", "sb441-n16.bsp", "earth_1962_240827_2124_combined.bpc", "earth_latest_high_prec.bpc"] {
        k.load(dir.join(f)).unwrap();
    }
    k
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let text = std::fs::read_to_string(&args[1]).unwrap();
    let recompute = args.get(3).map(|s| s == "spacerocks").unwrap_or(false);
    let k = kernel();

    let mut ids: Vec<String> = Vec::new();
    let mut groups: Vec<(Astrometry, Vec<String>)> = Vec::new();
    for line in text.lines().skip(1) {
        let c: Vec<&str> = line.split(',').collect();
        if ids.last().map(|s| s != c[0]).unwrap_or(true) {
            ids.push(c[0].to_string());
            groups.push((Astrometry::default(), Vec::new()));
        }
        let (a, stn) = groups.last_mut().unwrap();
        let v = |i: usize| c[i].parse::<f64>().unwrap();
        a.epoch.push(v(1));
        a.ra.push(v(2));
        a.dec.push(v(3));
        a.observer.push([v(4), v(5), v(6)]);
        a.sigma_ra.push(spacerocks::orbfit::DEFAULT_SIGMA);
        a.sigma_dec.push(spacerocks::orbfit::DEFAULT_SIGMA);
        stn.push(c[7].to_string());
    }
    if recompute {
        for (a, stn) in groups.iter_mut() {
            a.observer = observer_positions(stn, &a.epoch, &k).unwrap();
        }
    }

    let opts = FitOptions::default();
    let t0 = Instant::now();
    let fits: Vec<_> = groups
        .par_iter()
        .map(|(a, _)| {
            let t = Instant::now();
            let f = determine_orbit(a, None, &k, &opts);
            (f, t.elapsed().as_secs_f64())
        })
        .collect();
    eprintln!("{} objects in {:.2} s", fits.len(), t0.elapsed().as_secs_f64());

    let mut out = std::fs::File::create(&args[2]).unwrap();
    write!(out, "id,flag,csq,ndof,niter,epoch_tdb,x,y,z,xdot,ydot,zdot").unwrap();
    for i in 0..6 {
        for j in 0..6 {
            write!(out, ",cov_{}_{}", i, j).unwrap();
        }
    }
    writeln!(out, ",seconds").unwrap();
    for (id, (f, secs)) in ids.iter().zip(&fits) {
        write!(out, "{},{},{:?},{},{},{:?}", id, f.flag.code(), f.chi2, f.ndof, f.niter, f.epoch).unwrap();
        for v in f.state {
            write!(out, ",{:?}", v).unwrap();
        }
        for row in f.state_covariance() {
            for v in row {
                write!(out, ",{:?}", v).unwrap();
            }
        }
        writeln!(out, ",{:.3}", secs).unwrap();
    }
}
