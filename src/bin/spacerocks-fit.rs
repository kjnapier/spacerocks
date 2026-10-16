//! `spacerocks-fit`: fit an orbit to an object's MPC observations in the live dashboard.
//! See `spacerocks::orbfit::cli` (or run `spacerocks-fit --help`).
use spacerocks::orbfit::cli;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let open = |url: &str| cli::open_in_browser(url).map_err(|e| e.to_string());
    if let Err(e) = cli::run(&args, &open, None) {
        eprintln!("spacerocks-fit: {e}");
        std::process::exit(1);
    }
}
