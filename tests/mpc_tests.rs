//! Reading Minor Planet Center observations. Nothing here touches the network: the records and
//! API responses are written out by hand. Placing the stations needs a planetary ephemeris in
//! `SPACEROCKS_KERNELS`, with
//! Earth's orientation, and is skipped otherwise.

use std::path::PathBuf;

use serde_json::json;

use spacerocks::orbfit::mpc::{obs80_from_response, parse_obs80};
use spacerocks::orbfit::{mpc_astrometry, veres_sigma, ARCSEC};
use spacerocks::SpiceKernel;

/// Two usable records (G96 with Gaia DR2, 568 with Gaia DR3 under program code 2), then a B1950
/// one, a satellite one, an unknown station, a roving observer and a short line, which are
/// skipped, and a provisional designation from F51 without a catalog.
const OBS80: &str = "\
00433         C2019 03 15.47285 12 34 56.789-05 06 07.89         18.5 VV     G96
00433        2C2020 11 01.12345 00 00 01.500+89 59 59.99         18.5 VW     568
00433         A1940 01 01.5     12 34 56.789-05 06 07.89         18.5 VV     G96
00433         S2019 03 15.47285 12 34 56.789-05 06 07.89         18.5 VV     C51
00433         C2019 03 15.47285 12 34 56.789-05 06 07.89         18.5 VV     ZZZ
00433         C2019 03 15.47285 12 34 56.789-05 06 07.89         18.5 VV     247
00433         C2019 03 15.47285 12 34 56.789
     K20A00A   1999 12 31.99999 23 59 59.999-89 00 00.0          18.5 V      F51
";

#[test]
fn obs80_records() {
    let obs = parse_obs80(OBS80);
    assert_eq!(obs.len(), 3);
    let o = &obs[0];
    // 2019 March 15.47285 UTC.
    assert!((o.jd_utc - 2458557.97285).abs() < 1e-9, "{}", o.jd_utc);
    assert!((o.ra.to_degrees() - 15.0 * (12.0 + 34.0 / 60.0 + 56.789 / 3600.0)).abs() < 1e-12);
    assert!((o.dec.to_degrees() + (5.0 + 6.0 / 60.0 + 7.89 / 3600.0)).abs() < 1e-12);
    assert_eq!((o.station.as_str(), o.catalog.as_str(), o.program.as_str()), ("G96", "V", ""));
    assert_eq!((obs[1].station.as_str(), obs[1].catalog.as_str(), obs[1].program.as_str()), ("568", "W", "2"));
    assert!((obs[1].dec.to_degrees() - (90.0 - 0.01 / 3600.0)).abs() < 1e-12);
    // A blank note 2 is an ordinary record.
    assert_eq!((obs[2].station.as_str(), obs[2].catalog.as_str()), ("F51", ""));
    assert!((obs[2].dec.to_degrees() + 89.0).abs() < 1e-12);
    assert!((obs[2].jd_utc - 2451544.49999).abs() < 1e-9);
    assert!(parse_obs80("-00 30 00.00").is_empty());
}

#[test]
fn negative_zero_declinations() {
    let line = "00433         C2019 03 15.47285 12 34 56.789-00 30 00.0          18.5 VV     G96";
    let obs = parse_obs80(line);
    assert!((obs[0].dec.to_degrees() + 0.5).abs() < 1e-12);
}

#[test]
fn api_responses() {
    let text = "00433         C2019 03 15.47285 12 34 56.789-05 06 07.89         18.5 VV     G96";
    assert_eq!(obs80_from_response(&json!([{ "OBS80": text }])).unwrap(), text);
    assert_eq!(obs80_from_response(&json!([{ "OBS80": { "data": text } }])).unwrap(), text);
    assert_eq!(obs80_from_response(&json!([{ "OBS80": [text, text] }])).unwrap(), format!("{text}\n{text}"));
    assert!(obs80_from_response(&json!([])).is_err());
    assert!(obs80_from_response(&json!([{ "ADES_DF": [] }])).is_err());
    assert!(obs80_from_response(&json!({ "error": "unknown designation" })).unwrap_err().contains("unknown designation"));
}

#[test]
fn stations_placed_and_weighted() {
    let Some(dir) = std::env::var("SPACEROCKS_KERNELS").ok().map(PathBuf::from) else { return };
    let mut k = SpiceKernel::new();
    let Some(planets) = ["de440.bsp", "de440s.bsp", "de421.bsp"].into_iter().find(|f| dir.join(f).exists()) else { return };
    k.load(dir.join(planets)).unwrap();
    for pck in ["earth.tpc", "pck00011.tpc"] {
        if dir.join(pck).exists() {
            k.load(dir.join(pck)).unwrap();
        }
    }
    let obs = parse_obs80(OBS80);
    let a = match mpc_astrometry(&obs, &k, 1000) {
        Ok(a) => a,
        // Ground stations need Earth's orientation (an Earth binary PCK, or a text PCK).
        Err(e) if e.contains("orientation") => return,
        Err(e) => panic!("{e}"),
    };
    assert_eq!(a.len(), 3);
    // In time order, TDB 64-69 s after UTC.
    assert!(a.epoch.windows(2).all(|w| w[0] < w[1]));
    assert!(((a.epoch[0] - obs[2].jd_utc) * 86400.0 - 64.2).abs() < 1.0);
    // Vereš et al. weights; F51 gets 0.2".
    assert!((a.sigma_ra[0] / ARCSEC - 0.2).abs() < 1e-12);
    assert!((a.sigma_ra[1] / ARCSEC - veres_sigma("G96", a.epoch[1], Some("V"), None)).abs() < 1e-12);
    // Ground stations sit about an Earth radius from the geocenter.
    let earth = k.state_au(399, 0, a.epoch[1]).unwrap();
    let d = ((0..3).map(|j| (a.observer[1][j] - earth[j]).powi(2)).sum::<f64>()).sqrt() * 1.495978707e8;
    assert!((6300.0..6400.0).contains(&d), "{d} km");
    // Thinning keeps the ends of the arc.
    let thin = mpc_astrometry(&obs, &k, 2).unwrap();
    assert_eq!(thin.epoch, vec![a.epoch[0], a.epoch[2]]);
}
