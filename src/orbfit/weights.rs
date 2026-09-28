//! Astrometric uncertainties: layup's rendition of Vereš et al. (2017).

/// MPC one-character star-catalog codes (Obs80 column 72) for the catalogs the model branches on,
/// as their ADES names. Other codes and names pass through unchanged.
fn catalog_name(catalog: &str) -> &str {
    match catalog {
        "o" => "USNOB1",
        "s" => "USNOB2",
        "q" => "UCAC4",
        "t" => "PPMXL",
        "U" => "Gaia1",
        "V" => "Gaia2",
        "W" => "Gaia3",
        "X" => "Gaia3E",
        other => other,
    }
}

const EARLY_703_JD: f64 = 2456658.5;
const EARLY_691_JD: f64 = 2452640.5;
const EARLY_644_JD: f64 = 2452883.5;

/// One-sigma astrometric uncertainty (arcseconds, per axis) of an optical detection, from the
/// station, the TDB Julian date, and optionally the star catalog (ADES `astCat` name or MPC
/// one-letter code) and MPC program code: layup's `astrometric_uncertainty_Veres2017`, after
/// Vereš et al. (2017, Icarus 296, 139), reproduced branch for branch.
///
/// An empty `catalog` or `program` counts as absent. Stations not in the table get 1" with a
/// catalog and 1.5" without.
pub fn veres_sigma(station: &str, jd_tdb: f64, catalog: Option<&str>, program: Option<&str>) -> f64 {
    let catalog = catalog.filter(|c| !c.is_empty()).map(catalog_name);
    let gaia = |c: &str| matches!(c, "Gaia1" | "Gaia2" | "Gaia3" | "Gaia3E");
    let program = program.unwrap_or("");
    match station {
        "703" => {
            if jd_tdb <= EARLY_703_JD {
                1.0
            } else {
                0.8
            }
        }
        "691" => {
            if jd_tdb <= EARLY_691_JD {
                0.6
            } else {
                0.5
            }
        }
        "644" => {
            if jd_tdb <= EARLY_644_JD {
                0.6
            } else {
                0.4
            }
        }
        "704" => 1.0,
        "G96" => 0.5,
        "F51" => 0.2,
        "G45" => 0.6,
        "699" => 0.8,
        "D29" => 0.75,
        "C51" => 1.0,
        "E12" => 0.75,
        "608" => 0.6,
        "J75" => 1.0,
        "645" => 0.3,
        "673" => 0.3,
        "689" => 0.5,
        "950" => 0.5,
        "H01" => 0.3,
        "J04" => 0.4,
        // layup lists W84 here and again in the 0.4" group below; the first match wins.
        "W84" => 0.5,
        "G83" if program == "2" => match catalog {
            Some("UCAC4" | "PPMXL") => 0.3,
            Some(c) if gaia(c) => 0.2,
            Some(_) => 0.3,
            None => 1.0,
        },
        "K92" | "K93" | "Q63" | "Q64" | "V37" | "W85" | "W86" | "W87" | "K91" | "E10" | "F65" => 0.4,
        "Y28" => match catalog {
            Some("PPMXL" | "Gaia1") => 0.3,
            _ => 1.5,
        },
        "568" => match catalog {
            Some("USNOB1" | "USNOB2") => 0.5,
            Some(c) if gaia(c) => 0.1,
            Some("PPMXL") => 0.2,
            _ => 1.5,
        },
        "T09" | "T12" | "T14" => match catalog {
            Some(c) if gaia(c) => 0.1,
            _ => 1.5,
        },
        // Micheli's program; other catalogs keep layup's starting value, 1.5".
        "309" if program == "&" => match catalog {
            Some("UCAC4" | "PPMXL") => 0.3,
            Some(c) if gaia(c) => 0.2,
            _ => 1.5,
        },
        _ => {
            if catalog.is_some() {
                1.0
            } else {
                1.5
            }
        }
    }
}
