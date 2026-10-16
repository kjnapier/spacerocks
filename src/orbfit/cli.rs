//! `spacerocks-fit`: open the orbit-fit dashboard on an object's MPC observations from the
//! command line. The Rust binary of that name and the Python console script both run [`run`].
//!
//! ```text
//! spacerocks-fit                 # a random numbered asteroid
//! spacerocks-fit 433             # Eros
//! spacerocks-fit "2020 AA" --no-asteroids
//! ```
use std::path::PathBuf;

use serde_json::Value;

use crate::orbfit::{fetch_astrometry, FitDashboard, FitOptions, ASTEROID_PERTURBERS, LOAD_MAX_DETECTIONS};
use crate::spice::SpiceKernel;

/// The command line help.
pub const USAGE: &str = "\
usage: spacerocks-fit [designation] [options]

Fit an orbit to an object's observations from the Minor Planet Center in the live
orbit-fit dashboard. With no designation (or --random), picks a random numbered asteroid.
Load and Random on the page fetch other objects; press Done to print the fit and exit.

options:
  --random              a random numbered asteroid (the default without a designation)
  --port N              serve the page on 127.0.0.1:N (default: a free port)
  --no-browser          print the page's address without opening it
  --no-robust           start with robust outlier rejection off (MPC data spanning
                        decades usually needs it on)
  --no-asteroids        leave out the 16 massive asteroids (no sb441-n16.bsp needed)
  --max-detections N    keep at most N detections, spread over the arc (default 1000)
  --kernels FILE        load kernels from a TOML kernel-set file instead of the defaults
  -h, --help            show this help";

/// What the command line asked for.
#[derive(Debug, Clone, PartialEq)]
pub struct CliArgs {
    /// The object, or `None` for a random one.
    pub designation: Option<String>,
    pub port: u16,
    pub browser: bool,
    pub robust: bool,
    pub asteroids: bool,
    pub max_detections: usize,
    pub kernels: Option<PathBuf>,
}

impl Default for CliArgs {
    fn default() -> Self {
        CliArgs { designation: None, port: 0, browser: true, robust: true, asteroids: true, max_detections: LOAD_MAX_DETECTIONS, kernels: None }
    }
}

/// Read the arguments (without the program name). `Ok(None)` means help was asked for.
pub fn parse_args(args: &[String]) -> Result<Option<CliArgs>, String> {
    let mut a = CliArgs::default();
    let mut words: Vec<&str> = Vec::new();
    let mut random = false;
    let mut it = args.iter().map(String::as_str);
    while let Some(arg) = it.next() {
        let (flag, inline) = match arg.split_once('=') {
            Some((f, v)) if f.starts_with("--") => (f, Some(v.to_string())),
            _ => (arg, None),
        };
        let mut value = |name: &str| inline.clone().or_else(|| it.next().map(str::to_string)).ok_or(format!("{name} needs a value"));
        match flag {
            "-h" | "--help" => return Ok(None),
            "--random" => random = true,
            "--no-browser" => a.browser = false,
            "--no-robust" => a.robust = false,
            "--no-asteroids" => a.asteroids = false,
            "--port" => a.port = value("--port")?.parse().map_err(|_| "--port needs a port number")?,
            "--max-detections" => {
                a.max_detections = value("--max-detections")?.parse().map_err(|_| "--max-detections needs a whole number")?;
                if a.max_detections < 3 {
                    return Err("--max-detections must be at least 3".into());
                }
            }
            "--kernels" => a.kernels = Some(PathBuf::from(value("--kernels")?)),
            f if f.starts_with('-') => return Err(format!("unknown option {f}")),
            _ => words.push(arg),
        }
    }
    let designation = words.join(" ");
    if random && !designation.is_empty() {
        return Err("give a designation or --random, not both".into());
    }
    a.designation = (!designation.trim().is_empty()).then(|| designation.trim().to_string());
    Ok(Some(a))
}

/// Open `url` in the default browser (`open`, `xdg-open` or `start`).
pub fn open_in_browser(url: &str) -> std::io::Result<()> {
    use std::process::{Command, Stdio};
    let mut command = if cfg!(target_os = "macos") {
        Command::new("open")
    } else if cfg!(windows) {
        let mut c = Command::new("cmd");
        c.args(["/C", "start", ""]);
        c
    } else {
        Command::new("xdg-open")
    };
    command.arg(url).stdin(Stdio::null()).stdout(Stdio::null()).stderr(Stdio::null()).spawn().map(|_| ())
}

/// Run `spacerocks-fit` with `args` (without the program name): load the kernels, fetch the
/// object, serve the dashboard, open it with `open`, and once Done is pressed print the fit.
/// `interrupted` is checked while waiting; when it returns true this returns at once.
pub fn run(args: &[String], open: &dyn Fn(&str) -> Result<(), String>, interrupted: Option<&dyn Fn() -> bool>) -> Result<(), String> {
    let Some(a) = parse_args(args)? else {
        println!("{USAGE}");
        return Ok(());
    };
    eprintln!("loading kernels...");
    let kernel = match &a.kernels {
        Some(path) => SpiceKernel::from_config_file(path),
        None => SpiceKernel::defaults(),
    }
    .map_err(|e| format!("couldn't load the kernels: {e}"))?;
    match &a.designation {
        Some(d) => eprintln!("fetching {d} from the Minor Planet Center..."),
        None => eprintln!("fetching a random numbered asteroid from the Minor Planet Center..."),
    }
    let (name, astrometry) = fetch_astrometry(a.designation.as_deref(), &kernel, a.max_detections)?;
    eprintln!("{name}: {} detections", astrometry.len());

    let mut options = FitOptions { robust: a.robust, ..FitOptions::default() };
    if !a.asteroids {
        options.forces.exclude = ASTEROID_PERTURBERS.to_vec();
    }
    let max = a.max_detections;
    let dashboard = FitDashboard::start(a.port, &name, astrometry, kernel, options, None).map_err(|e| format!("couldn't start the dashboard: {e}"))?;
    dashboard.set_loader(move |designation, kernel| fetch_astrometry(designation, kernel, max));
    println!("orbit-fit dashboard: {}", dashboard.url());
    if a.browser {
        if let Err(e) = open(dashboard.url()) {
            eprintln!("couldn't open a browser ({e}); open the address above");
        }
    }
    eprintln!("press Done on the page to print the fit and exit (Ctrl-C quits)");
    if dashboard.wait(interrupted).is_none() {
        return Ok(());
    }
    println!("{}", report(&dashboard.name(), &dashboard.fit_summary()));
    Ok(())
}

/// The fit as text: outcome, epoch, elements with their uncertainties, and the state.
pub fn report(name: &str, fit: &Value) -> String {
    let f = |v: &Value| v.as_f64().unwrap_or(f64::NAN);
    if fit.is_null() {
        return format!("{name}: no fit");
    }
    let mut out = vec![format!(
        "{name}: {} · chi2/dof {:.4} ({} dof) · rms {:.3}\" · {} iterations",
        fit["status"].as_str().unwrap_or("?"),
        f(&fit["chi2_nu"]),
        fit["ndof"],
        f(&fit["rms"][2]),
        fit["niter"]
    )];
    out.push(format!("epoch      JD {:.6} TDB", f(&fit["epoch"])));
    let el = &fit["elements"];
    if el.is_object() {
        let rows = [("a (AU)", "a", 0), ("e", "e", 1), ("i (deg)", "inc", 2), ("node (deg)", "node", 3), ("peri (deg)", "arg", 4), ("M (deg)", "mean_anomaly", 5), ("q (AU)", "q", 6)];
        for (label, key, k) in rows {
            out.push(format!("{label:<10} {:.10} ± {:.2e}", f(&el[key]), f(&el["sigma"][k])));
        }
        out.push("           heliocentric ecliptic J2000".into());
    }
    if let Some(s) = fit["state"].as_array() {
        let s: Vec<String> = s.iter().map(|v| format!("{:.15e}", f(v))).collect();
        out.push(format!("state      [{}] (barycentric J2000, AU and AU/day)", s.join(", ")));
    }
    out.join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(s: &[&str]) -> Vec<String> {
        s.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn arguments() {
        assert_eq!(parse_args(&[]).unwrap(), Some(CliArgs::default()));
        assert_eq!(parse_args(&args(&["--help"])).unwrap(), None);
        let a = parse_args(&args(&["2020", "AA", "--no-robust", "--port=8800", "--no-browser", "--max-detections", "50", "--no-asteroids"])).unwrap().unwrap();
        assert_eq!(a.designation.as_deref(), Some("2020 AA"));
        assert!(!a.robust && !a.browser && !a.asteroids);
        assert_eq!((a.port, a.max_detections), (8800, 50));
        assert_eq!(parse_args(&args(&["--random"])).unwrap().unwrap().designation, None);
        assert!(parse_args(&args(&["433", "--random"])).is_err());
        assert!(parse_args(&args(&["--port"])).is_err());
        assert!(parse_args(&args(&["--frobnicate"])).is_err());
    }
}
