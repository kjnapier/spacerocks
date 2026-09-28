use numpy::IntoPyArray;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use spacerocks::orbfit;

use crate::rockcollection::extract_times;

/// A column of optional strings: None, a string for every row, or a sequence (None, NaN and ""
/// entries count as absent).
pub(crate) fn optional_strings(obj: Option<&Bound<'_, PyAny>>, n: usize, what: &str) -> PyResult<Vec<Option<String>>> {
    let Some(obj) = obj.filter(|o| !o.is_none()) else { return Ok(vec![None; n]) };
    if let Ok(s) = obj.extract::<String>() {
        return Ok(vec![Some(s); n]);
    }
    let items: Vec<Option<String>> = obj
        .try_iter()
        .map_err(|_| PyValueError::new_err(format!("{} must be a string or a sequence of strings", what)))?
        .map(|x| x.map(|x| x.extract::<String>().ok()))
        .collect::<PyResult<_>>()?;
    if items.len() != n {
        return Err(PyValueError::new_err(format!("{} has {} entries for {} detections", what, items.len(), n)));
    }
    Ok(items)
}

/// One-sigma astrometric uncertainties in radians (per axis, on-sky), after Vereš et al. (2017)
/// as layup assigns them (`weight_data=True`): by station, date, and optionally star catalog
/// (ADES astCat name or MPC one-letter code) and MPC program code. Pass the result as both
/// `sigma_ra` and `sigma_dec` to `fit`.
///
/// station: an MPC code or a sequence of them. epoch: Times, or Julian dates in `timescale`.
/// catalog, program: None, a string, or one per detection (None/NaN/"" = not given).
#[pyfunction]
#[pyo3(signature = (station, epoch, catalog=None, program=None, timescale="utc"))]
pub fn veres_sigma<'py>(
    py: Python<'py>,
    station: &Bound<'py, PyAny>,
    epoch: &Bound<'py, PyAny>,
    catalog: Option<&Bound<'py, PyAny>>,
    program: Option<&Bound<'py, PyAny>>,
    timescale: &str,
) -> PyResult<Bound<'py, numpy::PyArray1<f64>>> {
    let times = extract_times(epoch, timescale)?;
    let n = times.len();
    let stations: Vec<String> = if let Ok(s) = station.extract::<String>() {
        vec![s; n]
    } else {
        station.extract().map_err(|_| PyValueError::new_err("station must be an MPC code or a sequence of them"))?
    };
    if stations.len() != n {
        return Err(PyValueError::new_err(format!("{} stations for {} epochs", stations.len(), n)));
    }
    let catalogs = optional_strings(catalog, n, "catalog")?;
    let programs = optional_strings(program, n, "program")?;
    let arcsec = std::f64::consts::PI / (180.0 * 3600.0);
    let sigma: Vec<f64> = (0..n)
        .map(|i| orbfit::veres_sigma(&stations[i], times[i].tdb().jd(), catalogs[i].as_deref(), programs[i].as_deref()) * arcsec)
        .collect();
    Ok(sigma.into_pyarray(py))
}

static TABLES: std::sync::Mutex<Vec<(Option<std::path::PathBuf>, std::sync::Arc<orbfit::BiasTable>)>> = std::sync::Mutex::new(Vec::new());

/// A bias table, loaded once per process.
fn bias_table(path: Option<std::path::PathBuf>, download: bool) -> PyResult<std::sync::Arc<orbfit::BiasTable>> {
    let mut cache = TABLES.lock().unwrap();
    if let Some((_, t)) = cache.iter().find(|(p, _)| *p == path) {
        return Ok(t.clone());
    }
    let table = match &path {
        Some(p) => orbfit::BiasTable::load(p),
        None => orbfit::BiasTable::load_default(download),
    }
    .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let table = std::sync::Arc::new(table);
    cache.push((path, table.clone()));
    Ok(table)
}

/// Star-catalog debiasing (Eggl et al. 2020), as layup applies it with `debias_data=True`.
/// Returns the debiased (ra, dec), in radians.
///
/// ra, dec: radians. epoch: Times, or Julian dates in `timescale`. catalog: the star catalog of
/// each detection (ADES astCat name, e.g. "UCAC4", or MPC one-letter code, e.g. "q"): a string,
/// or one per detection; None/NaN/"" and catalogs the table doesn't cover (Gaia DR2/DR3/EDR3,
/// UCAC-5, ...) leave a detection unchanged. table: the path of JPL's `bias.dat`; by default the
/// one in `$SPACEROCKS_DEBIAS_DIR` or `~/.spacerocks/debias`, downloaded from JPL (~ 170 MB
/// compressed) if missing and `download` is true. The first load writes a binary copy
/// (`bias.bin`, ~330 MB) next to it; later loads memory-map that.
#[pyfunction]
#[pyo3(signature = (ra, dec, epoch, catalog, table=None, timescale="utc", download=true))]
#[allow(clippy::too_many_arguments)]
pub fn debias<'py>(
    py: Python<'py>,
    ra: Vec<f64>,
    dec: Vec<f64>,
    epoch: &Bound<'py, PyAny>,
    catalog: &Bound<'py, PyAny>,
    table: Option<std::path::PathBuf>,
    timescale: &str,
    download: bool,
) -> PyResult<(Bound<'py, numpy::PyArray1<f64>>, Bound<'py, numpy::PyArray1<f64>>)> {
    let times = extract_times(epoch, timescale)?;
    let n = ra.len();
    if dec.len() != n || times.len() != n {
        return Err(PyValueError::new_err(format!("ra, dec and epoch have lengths {}, {} and {}", n, dec.len(), times.len())));
    }
    let catalogs = optional_strings(Some(catalog), n, "catalog")?;
    let table = py.detach(|| bias_table(table, download))?;
    let jd: Vec<f64> = times.iter().map(|t| t.tdb().jd()).collect();
    let (mut ra_out, mut dec_out) = (ra.clone(), dec.clone());
    for i in 0..n {
        if let Some(c) = catalogs[i].as_deref() {
            let (r, d) = table.debias(ra[i], dec[i], jd[i], c);
            ra_out[i] = r;
            dec_out[i] = d;
        }
    }
    Ok((ra_out.into_pyarray(py), dec_out.into_pyarray(py)))
}

/// RA and Dec (radians) of an occulting object from ADES occultation astrometry: the occulted
/// star's position `ra_star`, `dec_star` and the object's offset from it, `delta_ra` (which
/// includes cos Dec) and `delta_dec`, all in radians (ADES gives the star in degrees and the
/// offsets in arcsec). Scalars or arrays. The offset is applied on the tangent plane at the star.
#[pyfunction]
pub fn occultation_radec<'py>(
    py: Python<'py>,
    ra_star: &Bound<'py, PyAny>,
    dec_star: &Bound<'py, PyAny>,
    delta_ra: &Bound<'py, PyAny>,
    delta_dec: &Bound<'py, PyAny>,
) -> PyResult<(Bound<'py, numpy::PyArray1<f64>>, Bound<'py, numpy::PyArray1<f64>>)> {
    let col = |o: &Bound<'py, PyAny>| -> PyResult<Vec<f64>> {
        if let Ok(v) = o.extract::<f64>() {
            return Ok(vec![v]);
        }
        o.extract::<Vec<f64>>()
    };
    let (a, d, da, dd) = (col(ra_star)?, col(dec_star)?, col(delta_ra)?, col(delta_dec)?);
    let n = a.len().max(d.len()).max(da.len()).max(dd.len());
    let get = |v: &Vec<f64>, i: usize, what: &str| -> PyResult<f64> {
        match v.len() {
            1 => Ok(v[0]),
            m if m == n => Ok(v[i]),
            m => Err(PyValueError::new_err(format!("{} has {} entries, expected 1 or {}", what, m, n))),
        }
    };
    let mut ra = Vec::with_capacity(n);
    let mut dec = Vec::with_capacity(n);
    for i in 0..n {
        let (r, s) = orbfit::occultation_radec(get(&a, i, "ra_star")?, get(&d, i, "dec_star")?, get(&da, i, "delta_ra")?, get(&dd, i, "delta_dec")?);
        ra.push(r);
        dec.push(s);
    }
    Ok((ra.into_pyarray(py), dec.into_pyarray(py)))
}
