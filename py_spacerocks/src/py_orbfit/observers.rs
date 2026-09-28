use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;

use spacerocks::orbfit;

use crate::py_orbfit::weights::optional_strings;
use crate::rockcollection::extract_times;
use crate::PySpiceKernel;

/// Rows of an optional (n, 3) array, NaN rows (or no array) meaning "not given".
fn optional_rows(obj: Option<PyReadonlyArray2<'_, f64>>, n: usize, what: &str) -> PyResult<Vec<Option<[f64; 3]>>> {
    let Some(arr) = obj else { return Ok(vec![None; n]) };
    let a = arr.as_array();
    if a.shape() != [n, 3] {
        return Err(PyValueError::new_err(format!("{} has shape {:?}, expected ({}, 3)", what, a.shape(), n)));
    }
    Ok((0..n)
        .map(|i| {
            let r = [a[[i, 0]], a[[i, 1]], a[[i, 2]]];
            if r.iter().all(|v| v.is_finite()) { Some(r) } else { None }
        })
        .collect())
}

/// Barycentric J2000 states (AU, AU/day), shape (n, 6), of the observers of detections, for
/// `fit(observer=...)`. A detection with a position (`pos`, as ADES pos1..pos3) is placed from
/// it, as layup does: `sys` "ICRF_KM" or "ICRF_AU" is a geocentric ICRF position (km or AU), with
/// the optional velocity `vel` (km/s or AU/day; otherwise the Earth's), and "WGS84" is east
/// longitude and geodetic latitude (degrees) and height (m), e.g. a roving observer (247), which
/// turns with the Earth. `ctr` must be 399 (the default). The other detections are placed from
/// their MPC station codes: fixed stations from their parallax constants, spacecraft from a
/// loaded SPK kernel or JPL Horizons.
///
/// station: MPC code(s). epoch: Times, or Julian dates in `timescale`. sys: None, a string, or
/// one per detection. ctr: None, an int, or one per detection. pos, vel: (n, 3) arrays with NaN
/// rows where not given.
#[pyfunction]
#[pyo3(signature = (station, epoch, kernel, sys=None, ctr=None, pos=None, vel=None, timescale="utc"))]
#[allow(clippy::too_many_arguments)]
pub fn observers<'py>(
    py: Python<'py>,
    station: &Bound<'py, PyAny>,
    epoch: &Bound<'py, PyAny>,
    kernel: PyRef<PySpiceKernel>,
    sys: Option<&Bound<'py, PyAny>>,
    ctr: Option<&Bound<'py, PyAny>>,
    pos: Option<PyReadonlyArray2<'py, f64>>,
    vel: Option<PyReadonlyArray2<'py, f64>>,
    timescale: &str,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
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
    let systems = optional_strings(sys, n, "sys")?;
    let centers: Vec<i64> = match ctr.filter(|c| !c.is_none()) {
        None => vec![399; n],
        Some(c) => match c.extract::<i64>() {
            Ok(v) => vec![v; n],
            Err(_) => c
                .try_iter()?
                .map(|x| x.and_then(|x| Ok(x.extract::<f64>().ok().filter(|v| v.is_finite()).map(|v| v as i64).unwrap_or(399))))
                .collect::<PyResult<_>>()?,
        },
    };
    if centers.len() != n {
        return Err(PyValueError::new_err(format!("ctr has {} entries for {} detections", centers.len(), n)));
    }
    let positions = optional_rows(pos, n, "pos")?;
    let velocities = optional_rows(vel, n, "vel")?;
    let jd: Vec<f64> = times.iter().map(|t| t.tdb().jd()).collect();
    let k = &kernel.inner;

    let given: Vec<usize> = (0..n).filter(|&i| positions[i].is_some()).collect();
    let coded: Vec<usize> = (0..n).filter(|&i| positions[i].is_none()).collect();
    let out = py.detach(|| -> Result<Vec<[f64; 6]>, String> {
        let mut out = vec![[f64::NAN; 6]; n];
        let from_pos: Vec<[f64; 6]> = given
            .par_iter()
            .map(|&i| {
                let sys = systems[i].as_deref().ok_or_else(|| format!("detection {} has a position but no sys", i))?;
                orbfit::ades_observer_state(sys, centers[i], positions[i].unwrap(), velocities[i], jd[i], k).map_err(|e| format!("detection {}: {}", i, e))
            })
            .collect::<Result<_, String>>()?;
        let codes: Vec<String> = coded.iter().map(|&i| stations[i].clone()).collect();
        let epochs: Vec<f64> = coded.iter().map(|&i| jd[i]).collect();
        let from_code = orbfit::observer_states(&codes, &epochs, k)?;
        for (s, &i) in from_pos.iter().zip(&given) {
            out[i] = *s;
        }
        for (s, &i) in from_code.iter().zip(&coded) {
            out[i] = *s;
        }
        Ok(out)
    })
    .map_err(PyValueError::new_err)?;
    Ok(numpy::ndarray::Array2::from_shape_fn((n, 6), |(i, j)| out[i][j]).into_pyarray(py))
}
