use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use spacerocks::spice::{et_from_jd, jd_from_et, km_to_au_state, KernelConfig, KernelKind, SpiceKernel};
use pyo3::types::{PyDict, PyType};

use crate::py_time::time::PyTime;

fn err<E: std::fmt::Display>(e: E) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// A self-contained collection of loaded SPICE kernels.
///
/// Every SpiceKernel is independent: loading into one does not affect another. `copy()`
/// returns a new kernel that shares the already-loaded (memory-mapped) files, so it is cheap.
#[pyclass]
#[pyo3(name = "SpiceKernel")]
pub struct PySpiceKernel {
    pub inner: SpiceKernel,
}

impl PySpiceKernel {
    fn resolve_body(&self, body: &Bound<'_, PyAny>) -> PyResult<i32> {
        if let Ok(id) = body.extract::<i32>() {
            return Ok(id);
        }
        let name: String = body.extract()?;
        self.inner
            .body_id(&name)
            .ok_or_else(|| PyValueError::new_err(format!("unknown body '{}'", name)))
    }
}

#[pymethods]
impl PySpiceKernel {
    #[new]
    fn new() -> Self {
        PySpiceKernel { inner: SpiceKernel::new() }
    }

    /// A kernel with the standard spacerocks kernel set loaded.
    ///
    /// Kernels are looked up in the cache directory (`~/.spacerocks/spice`, or the
    /// `SPACEROCKS_SPICE_DIR` environment variable) and downloaded from NAIF/JPL if missing:
    /// `latest_leapseconds.tls`, `de440s.bsp`, the newest `earth_1962_*_combined.bpc`,
    /// `gm_de440.tpc`, and `sb441-n16.bsp` (asteroid perturbers, ~650 MB on first use).
    ///
    /// `download=False` uses only local files. `update=True` checks for a newer Earth
    /// orientation file even if one is cached.
    #[classmethod]
    #[pyo3(signature = (download = true, update = false))]
    fn defaults(_cls: &Bound<'_, PyType>, py: Python<'_>, download: bool, update: bool) -> PyResult<Self> {
        let mut cfg = KernelConfig::default_with_download(download);
        cfg.check_for_updates = update;
        let k = py.detach(|| SpiceKernel::from_config(&cfg)).map_err(err)?;
        Ok(PySpiceKernel { inner: k })
    }

    /// A kernel built from a TOML kernel-set file (see the spacerocks docs for the format).
    /// `download` overrides the file's `auto_download` setting when given.
    #[classmethod]
    #[pyo3(signature = (path, download = None))]
    fn from_config(_cls: &Bound<'_, PyType>, py: Python<'_>, path: &str, download: Option<bool>) -> PyResult<Self> {
        let mut cfg = KernelConfig::from_file(path).map_err(err)?;
        if let Some(d) = download {
            cfg.auto_download = d;
        }
        let k = py.detach(|| SpiceKernel::from_config(&cfg)).map_err(err)?;
        Ok(PySpiceKernel { inner: k })
    }

    /// Load the standard kernel set into this kernel (see `defaults`).
    #[pyo3(signature = (download = true, update = false))]
    fn load_defaults(&mut self, py: Python<'_>, download: bool, update: bool) -> PyResult<()> {
        let mut cfg = KernelConfig::default_with_download(download);
        cfg.check_for_updates = update;
        let inner = &mut self.inner;
        py.detach(|| inner.load_config(&cfg)).map_err(err)
    }

    /// The directory where downloaded kernels are cached.
    #[staticmethod]
    fn cache_dir() -> String {
        spacerocks::spice::config::default_download_dir().display().to_string()
    }

    /// Load any kernel: SPK, binary PCK, text kernel (LSK, text PCK, FK, ...) or meta-kernel.
    fn load(&mut self, path: &str) -> PyResult<()> {
        self.inner.load(path).map_err(err)
    }

    fn load_spk(&mut self, path: &str) -> PyResult<()> {
        self.inner.load_spk(path).map_err(err)
    }

    fn load_bpc(&mut self, path: &str) -> PyResult<()> {
        self.inner.load_pck(path).map_err(err)
    }

    fn load_pck(&mut self, path: &str) -> PyResult<()> {
        self.inner.load_pck(path).map_err(err)
    }

    /// Unload a kernel (a meta-kernel also unloads what it loaded). Returns True if found.
    fn unload(&mut self, path: &str) -> bool {
        self.inner.unload(path)
    }

    /// Unload everything.
    fn clear(&mut self) {
        self.inner.clear();
    }

    /// A new kernel sharing this kernel's loaded files.
    fn copy(&self) -> Self {
        PySpiceKernel { inner: self.inner.clone() }
    }

    fn __copy__(&self) -> Self {
        self.copy()
    }

    /// List of (path, kind) for loaded files, lowest priority first.
    #[getter]
    fn loaded_kernels(&self) -> Vec<(String, String)> {
        self.inner
            .loaded_kernels()
            .into_iter()
            .map(|k| (k.path, k.kind.to_string()))
            .collect()
    }

    /// NAIF IDs of all bodies with SPK data.
    #[getter]
    fn bodies(&self) -> Vec<i32> {
        self.inner.spk_bodies()
    }

    fn body_id(&self, name: &str) -> Option<i32> {
        self.inner.body_id(name)
    }

    fn body_name(&self, id: i32) -> Option<String> {
        self.inner.body_name(id)
    }

    /// SPK coverage of a body as a list of (start, end) TDB Julian dates.
    fn coverage(&self, body: &Bound<'_, PyAny>) -> PyResult<Vec<(f64, f64)>> {
        let id = self.resolve_body(body)?;
        Ok(self
            .inner
            .spk_coverage(id)
            .into_iter()
            .map(|iv| (jd_from_et(iv.start), jd_from_et(iv.end)))
            .collect())
    }

    /// Geometric state of `target` relative to `observer` (names or NAIF IDs) at `epoch`.
    ///
    /// Returns a length-6 array: position and velocity in `frame`, in AU and AU/day
    /// (`units="au"`, the default) or km and km/s (`units="km"`).
    #[pyo3(signature = (target, observer, epoch, frame = "J2000", units = "au"))]
    fn state<'py>(
        &self,
        py: Python<'py>,
        target: &Bound<'_, PyAny>,
        observer: &Bound<'_, PyAny>,
        epoch: PyRef<PyTime>,
        frame: &str,
        units: &str,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let t = self.resolve_body(target)?;
        let o = self.resolve_body(observer)?;
        let et = et_from_jd(epoch.inner.tdb().jd());
        let s = self.inner.state(t, o, frame, et).map_err(err)?;
        let s = match units.to_lowercase().as_str() {
            "au" => km_to_au_state(&s),
            "km" => s,
            other => return Err(PyValueError::new_err(format!("units must be 'au' or 'km', not '{}'", other))),
        };
        Ok(s.to_vec().into_pyarray(py))
    }

    /// 3x3 rotation matrix taking vectors from frame `from_frame` to `to_frame` at `epoch`.
    fn pxform<'py>(&self, py: Python<'py>, from_frame: &str, to_frame: &str, epoch: PyRef<PyTime>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let et = et_from_jd(epoch.inner.tdb().jd());
        let x = self.inner.sxform(from_frame, to_frame, et).map_err(err)?;
        let arr = numpy::ndarray::Array2::from_shape_fn((3, 3), |(i, j)| x.rotation[i][j]);
        Ok(arr.into_pyarray(py))
    }

    /// 6x6 state transformation matrix from `from_frame` to `to_frame` at `epoch`
    /// (derivative block per second, as in SPICE).
    fn sxform<'py>(&self, py: Python<'py>, from_frame: &str, to_frame: &str, epoch: PyRef<PyTime>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let et = et_from_jd(epoch.inner.tdb().jd());
        let x = self.inner.sxform(from_frame, to_frame, et).map_err(err)?;
        let arr = numpy::ndarray::Array2::from_shape_fn((6, 6), |(i, j)| match (i < 3, j < 3) {
            (true, true) => x.rotation[i][j],
            (false, false) => x.rotation[i - 3][j - 3],
            (false, true) => x.rate[i - 3][j],
            (true, false) => 0.0,
        });
        Ok(arr.into_pyarray(py))
    }

    fn __repr__(&self) -> String {
        super::display::text(&self.inner)
    }

    /// Rich notebook display: loaded files by priority, a coverage timeline per body and
    /// frame, and a segment table.
    fn _repr_html_(&self) -> String {
        super::display::html(&self.inner)
    }

    /// What each loaded file contains, lowest priority first, as a list of dicts with keys
    /// `path`, `kind`, `size_bytes`, `description`, `segments`, `variables`, `children`.
    /// Each segment group has `body`, `name`, `center`, `frame`, `type`, `n_segments` and
    /// `coverage` (list of (start, end) TDB Julian dates). For binary PCK files `body` is
    /// the frame class ID and `frame` the inertial frame the orientation is relative to.
    fn summary<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let k = &self.inner;
        k.summary()
            .iter()
            .map(|f| {
                let d = PyDict::new(py);
                d.set_item("path", &f.path)?;
                d.set_item("kind", f.kind.to_string())?;
                d.set_item("size_bytes", f.size_bytes)?;
                d.set_item("description", super::display::describe(k, f))?;
                let segs = f
                    .groups
                    .iter()
                    .map(|g| {
                        let s = PyDict::new(py);
                        s.set_item("body", g.body)?;
                        let name = if f.kind == KernelKind::Pck {
                            Some(super::display::pck_frame_label(k, g.body))
                        } else {
                            k.body_name(g.body)
                        };
                        s.set_item("name", name)?;
                        s.set_item("center", g.center)?;
                        s.set_item("frame", k.frame_info(g.frame).map(|x| x.name).unwrap_or_else(|| g.frame.to_string()))?;
                        s.set_item("type", g.data_type)?;
                        s.set_item("n_segments", g.n_segments)?;
                        let cov: Vec<(f64, f64)> = g.coverage.iter().map(|iv| (jd_from_et(iv.start), jd_from_et(iv.end))).collect();
                        s.set_item("coverage", cov)?;
                        Ok(s)
                    })
                    .collect::<PyResult<Vec<_>>>()?;
                d.set_item("segments", segs)?;
                d.set_item("variables", &f.variables)?;
                d.set_item("children", &f.children)?;
                Ok(d)
            })
            .collect()
    }

    fn __len__(&self) -> usize {
        self.inner.loaded_kernels().len()
    }
}
