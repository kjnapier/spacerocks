use pyo3::prelude::*;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::types::{PyDict, IntoPyDict, PyList, PyTuple};
use std::path::{Path, PathBuf};
use serde::{Deserialize, Serialize};
use flate2::read::GzDecoder;
use std::fs;
use std::io::Read;
use arrow::array::{Float64Array, StringArray, Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use rayon::prelude::*;

use spacerocks::spacerock::SpaceRock;
use spacerocks::Time;
use spacerocks::transforms::calc_true_anomaly_from_mean_anomaly;
use spacerocks::Observatory;
// use crate::py_time::time::PyTime;
use crate::py_time::time::PyTime;
use crate::py_observing::observer::PyObserver;
use crate::py_observing::observation::PyObservation;

use serde_json;

use crate::RockCollection;
use anyhow::Error;
use spacerocks::constants::MPC_URL;
use dirs::home_dir;
use std::time::Duration;
use reqwest::ClientBuilder;
use spacerocks::constants::MPC_API;
use spacerocks::observing::{Observer, Observation};


#[derive(Debug, Clone)]
pub enum StorageFormat {
    None,
    JsonGz,
    Feather,
}

#[derive(Debug, Clone)]
pub enum OutputFormat {
    DataFrame,
    RockCollection,
}


#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MPCData {
    pub H: Option<f64>, 
    pub G: Option<f64>,
    pub Epoch: f64, 
    pub M: f64,
    pub Peri: f64,
    pub Node: f64,
    pub i: f64,
    pub e: f64,
    pub a: f64,
    pub Principal_desig: String,
    #[serde(rename = "Orbit_type")]
    pub orbit_type: Option<String>,
}

#[pyclass]
pub struct MPCHandler {
    path: PathBuf,
}

#[pymethods]
impl MPCHandler {
    #[new]
    #[pyo3(signature = (path=None))]
    pub fn new(path: Option<PathBuf>) -> Self {
        let default_path = home_dir()
            .unwrap_or_default()
            .join(".spacerocks")
            .join("mpc");
        
        let final_path = path.unwrap_or(default_path);
        

        MPCHandler {
            path: final_path,
        }
    }

    /// Fetch and optionally process MPC data
    #[pyo3(signature = (catalog, orbit_type=None, storage_format=Some("feather".to_string()), output_format="dataframe"))]
    pub fn fetch_data(&self, catalog: String, orbit_type: Option<String>, storage_format: Option<String>, output_format: &str) -> PyResult<PyObject> {

        // Get the storage format
        let storage_type = match storage_format.as_deref() {
            Some("feather") => StorageFormat::Feather,
            Some("json") => StorageFormat::JsonGz,
            None => StorageFormat::None,
            Some(other) => return Err(PyErr::new::<PyValueError, _>(
                format!("Invalid storage format: {}. Must be 'feather', 'json', or None", other)
            )),
        };


        // Get output format
        let out_format = match output_format {
            "dataframe" => OutputFormat::DataFrame,
            "rocks" => OutputFormat::RockCollection,
            _ => return Err(PyErr::new::<PyValueError, _>(
                "Output format must be 'dataframe' or 'rocks'"
            )),
        };

        // Get the data
        let mut data = self.get_data(&catalog, &storage_type)?;

        // Filter by orbit_type if specified
        if let Some(orbit_filter) = orbit_type {
        //     data.retain(|d| {
        //         d.orbit_type.as_ref().map_or(false, |ot| ot == &orbit_filter)
        //     });
        // }
                let filtered: Vec<_> = data.par_iter()
                .filter(|d| d.orbit_type.as_ref().map_or(false, |ot| ot == &orbit_filter))
                .cloned()
                .collect();
            data = filtered;
        }

        // Convert to requested output format
        match out_format {
            OutputFormat::DataFrame => self.to_dataframe(data),
            OutputFormat::RockCollection => self.to_rock_collection(data),
        }
    }

    #[staticmethod]
    pub fn create_rock_collection(mpc_path: PathBuf, catalog: String, download_data: bool, orbit_type: Option<String>,) -> PyResult<RockCollection> {
        let handler = MPCHandler::new(Some(mpc_path));

        // Create the necessary directories
        fs::create_dir_all(&handler.path)?;

        // Only download if requested or if data does not exist
        if download_data {
            println!("Downloading new MPC data...");
            handler.fetch_data(
                catalog.clone(),
                orbit_type.clone(),
                Some("feather".to_string()),
                "rocks"
            )?;
        }

        // Convert PyObject to RockCollection
        let py_obj = handler.fetch_data(
            catalog,
            orbit_type,
            Some("feather".to_string()),
            "rocks"
        )?;

        Python::with_gil(|py| {
            let rocks_collection = py_obj.extract::<RockCollection>(py)?;
            Ok(rocks_collection)
        })
    }

    #[pyo3(signature = (designation, formats=vec!["ADES_DF".to_string()], ades_version=None))]
    pub fn get_detections(&self, designation: &str, formats: Vec<String>, ades_version: Option<String>) -> PyResult<PyObject> {
        Python::with_gil(|py| {
            let requests = PyModule::import(py, "requests")?;
            
            let request_data = PyDict::new(py);
            request_data.set_item("desigs", vec![designation.to_string()])?;
            request_data.set_item("output_format", formats)?;
            if let Some(version) = ades_version {
                request_data.set_item("ades_version", version)?;
            }

            let kwargs = PyDict::new(py);
            kwargs.set_item("json", request_data)?;

            let get_func = requests.getattr("get")?;
            let response = get_func.call((MPC_API,), Some(&kwargs))?;
            
            let json_response = response.call_method0("json")?;

            // Extract ADES data and create observations
            let ades_data = json_response.get_item(0)?.get_item("ADES_DF")?;
            let mut observations = Vec::new();

            // Iterate over each row in the ADES data
            for row in ades_data.try_iter()? {
                let row = row?;

                // Time conversion
                let time = Time::from_isot(&row.get_item("obstime")?.extract::<String>()?)
                    .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
                
                // Observatory setup
                let observatory = Observatory::from_obscode(row.get_item("stn")?.extract()?);
                let observer = observatory.expect("Failed at making observatory").at(&time, "J2000", "SSB");

                // Convert ra and dec to radians
                let ra = row.get_item("ra")?
                    .extract::<String>()?
                    .parse::<f64>()
                    .map_err(|_| pyo3::exceptions::PyValueError::new_err("Invalid RA format"))?
                    .to_radians();

                let dec = row.get_item("dec")?
                    .extract::<String>()?
                    .parse::<f64>()
                    .map_err(|_| pyo3::exceptions::PyValueError::new_err("Invalid Dec format"))?
                    .to_radians();

                let mag = match row.get_item("mag") {
                    Ok(m) => match m.extract::<String>() {
                        Ok(val) => val.parse::<f64>().ok(), // Convert string to f64
                        Err(_) => m.extract::<f64>().ok(),  // If already f64, use it directly
                    },
                    Err(_) => None,
                };


                let observation = PyObservation { 
                    inner: Observation::from_astrometry(
                        time,
                        ra, 
                        dec,
                        mag,
                        observer.expect("Couldn't make Observer object"),
                    )
                };

                observations.push(Py::new(py, observation)?);
            }

            let py_list = PyList::new(py, observations)?;

            let json_response: PyObject = json_response.to_object(py);
            let elements: &[PyObject] = &[py_list.to_object(py), json_response.to_object(py)];
            let tuple = PyTuple::new(py, elements);


            Ok(tuple?.to_object(py))
        })
    }

}

impl MPCHandler {

    fn get_data(&self, catalog: &str, storage_type: &StorageFormat) -> PyResult<Vec<MPCData>> {
        // First check if we already have the data in any format
        let feather_path = self.path.join(format!("{}.feather", catalog));
        let json_path = self.path.join(format!("{}.json.gz", catalog));
    
        // If data exists, use it regardless of specified storage_type
        if json_path.exists() {
            println!("Using existing json file: {}", json_path.display());
            let reader = GzDecoder::new(fs::File::open(&json_path)?);
            return serde_json::from_reader(reader)
                .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()));
        } else if feather_path.exists() {
            println!("Using existing feather file: {}", feather_path.display());
            return self.read_from_feather(&feather_path);
        }
    
        // If we get here, we need to download
        let url = format!(
            "{}/{}.json.gz",
            MPC_URL,
            catalog
        );
    
        let client = reqwest::blocking::ClientBuilder::new()
            .timeout(Duration::from_secs(300))  // 5 minute timeout
            .build()
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;
    
        // Download and process based on desired storage format
        match storage_type {
            StorageFormat::None => {
                println!("Downloading data from {}", url);
                let response = client
                    .get(&url)
                    .send()
                    .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;
                let reader = GzDecoder::new(response);
                serde_json::from_reader(reader)
                    .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))
            },
            StorageFormat::JsonGz => {
                println!("Downloading data from {}", url);
                println!("Saving to {}", json_path.display());
                let response = client
                    .get(&url)
                    .send()
                    .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;
                let content = response.bytes()
                    .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;
                fs::create_dir_all(&self.path)?;
                fs::write(&json_path, content)?;
                let reader = GzDecoder::new(fs::File::open(&json_path)?);
                serde_json::from_reader(reader)
                    .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))
            },
            StorageFormat::Feather => {
                println!("Downloading data from {}", url);
                println!("Converting and saving to {}", feather_path.display());
                let response = client
                    .get(&url)
                    .send()
                    .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;
                let reader = GzDecoder::new(response);
                let data: Vec<MPCData> = serde_json::from_reader(reader)
                    .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;
                
                fs::create_dir_all(&self.path)?;
                self.save_as_feather(&data, &feather_path)?;
                Ok(data)
            }
        }
    }

    // fn to_dataframe(&self, data: Vec<MPCData>) -> PyResult<PyObject> {
    //     Python::with_gil(|py| {
    //         let pd = py.import("pandas")?;
    //         let dict = PyDict::new(py);
    //         dict.set_item("H", data.iter().map(|d| d.H).collect::<Vec<_>>())?;
    //         dict.set_item("G", data.iter().map(|d| d.G).collect::<Vec<_>>())?;
    //         dict.set_item("Epoch", data.iter().map(|d| d.Epoch).collect::<Vec<_>>())?;
    //         dict.set_item("M", data.iter().map(|d| d.M).collect::<Vec<_>>())?;
    //         dict.set_item("Peri", data.iter().map(|d| d.Peri).collect::<Vec<_>>())?;
    //         dict.set_item("Node", data.iter().map(|d| d.Node).collect::<Vec<_>>())?;
    //         dict.set_item("i", data.iter().map(|d| d.i).collect::<Vec<_>>())?;
    //         dict.set_item("e", data.iter().map(|d| d.e).collect::<Vec<_>>())?;
    //         dict.set_item("a", data.iter().map(|d| d.a).collect::<Vec<_>>())?;
    //         dict.set_item("Principal_desig", data.iter().map(|d| &d.Principal_desig).collect::<Vec<_>>())?;
    //         dict.set_item("orbit_type", data.iter().map(|d| &d.orbit_type).collect::<Vec<_>>())?;
            
    //         let df = pd.call_method1("DataFrame", (dict,))?;
    //         Ok(df.into_py(py))
    //     })
    // }

    fn to_dataframe(&self, data: Vec<MPCData>) -> PyResult<PyObject> {
        Python::with_gil(|py| {
            let pd = py.import("pandas")?;
            let dict = PyDict::new(py);
            
            // First join: split into two main groups
            let (group1, group2) = rayon::join(
                || {
                    // Second join: split first group
                    rayon::join(
                        || {
                            // Process H, G, Epoch
                            (
                                data.par_iter().map(|d| d.H).collect::<Vec<_>>(),
                                data.par_iter().map(|d| d.G).collect::<Vec<_>>(),
                                data.par_iter().map(|d| d.Epoch).collect::<Vec<_>>()
                            )
                        },
                        || {
                            // Process M, Peri, Node
                            (
                                data.par_iter().map(|d| d.M).collect::<Vec<_>>(),
                                data.par_iter().map(|d| d.Peri).collect::<Vec<_>>(),
                                data.par_iter().map(|d| d.Node).collect::<Vec<_>>()
                            )
                        }
                    )
                },
                || {
                    // Third join: split second group
                    rayon::join(
                        || {
                            // Process i, e, a
                            (
                                data.par_iter().map(|d| d.i).collect::<Vec<_>>(),
                                data.par_iter().map(|d| d.e).collect::<Vec<_>>(),
                                data.par_iter().map(|d| d.a).collect::<Vec<_>>()
                            )
                        },
                        || {
                            // Process Principal_desig and orbit_type
                            (
                                data.par_iter().map(|d| &d.Principal_desig).collect::<Vec<_>>(),
                                data.par_iter().map(|d| &d.orbit_type).collect::<Vec<_>>()
                            )
                        }
                    )
                }
            );
    
            let ((h, g, epoch), (m, peri, node)) = group1;
            let ((inc, e, a), (desig, orbit)) = group2;
    
            // Set all items
            dict.set_item("H", h)?;
            dict.set_item("G", g)?;
            dict.set_item("Epoch", epoch)?;
            dict.set_item("M", m)?;
            dict.set_item("Peri", peri)?;
            dict.set_item("Node", node)?;
            dict.set_item("i", inc)?;
            dict.set_item("e", e)?;
            dict.set_item("a", a)?;
            dict.set_item("Principal_desig", desig)?;
            dict.set_item("orbit_type", orbit)?;
            
            let df = pd.call_method1("DataFrame", (dict,))?;
            Ok(df.into_py(py))
        })
    }

    fn to_rock_collection(&self, data: Vec<MPCData>) -> PyResult<PyObject> {
        let rocks: Vec<Result<SpaceRock, String>> = 
            data.par_iter().map(|d| {
                match SpaceRock::from_kepler(
                    &d.Principal_desig,
                    d.a * (1.0 - d.e), // q, not a
                    d.e,
                    d.i.to_radians(),
                    d.Peri.to_radians(),
                    d.Node.to_radians(),
                    calc_true_anomaly_from_mean_anomaly(d.e, d.M.to_radians()).expect("Error calculating true anomaly"),
                    Time::new(d.Epoch, "utc", "jd").map_err(|e| e.to_string())?,
                    "ECLIPJ2000",
                    "SSB",
                ) {
                    Ok(mut rock) => {
                        if let Some(h) = d.H {
                            rock.set_absolute_magnitude(h);
                        }
                        if let Some(g) = d.G {
                            rock.set_gslope(g);
                        }
                        Ok(rock)
                    },
                    Err(e) => Err(e.to_string())
                }
            }).collect();
    
        // Process results
        let rocks: Result<Vec<SpaceRock>, String> = rocks.into_iter()
            .collect();
    
        match rocks {
            Ok(rocks) => Python::with_gil(|py| {
                #[allow(deprecated)]
                Ok(RockCollection { rocks }.into_py(py))
            }),
            Err(e) => Err(PyErr::new::<PyRuntimeError, _>(e))
        }
    }

    fn save_as_feather(&self, data: &[MPCData], path: &Path) -> PyResult<()> {
        let schema = Schema::new(vec![
            Field::new("H", DataType::Float64, true),
            Field::new("G", DataType::Float64, true),
            Field::new("Epoch", DataType::Float64, true),  // Changed to nullable
            Field::new("M", DataType::Float64, true),      // Changed to nullable
            Field::new("Peri", DataType::Float64, true),   // Changed to nullable
            Field::new("Node", DataType::Float64, true),   // Changed to nullable
            Field::new("i", DataType::Float64, true),      // Changed to nullable
            Field::new("e", DataType::Float64, true),      // Changed to nullable
            Field::new("a", DataType::Float64, true),      // Changed to nullable
            Field::new("Principal_desig", DataType::Utf8, true), // Changed to nullable
            Field::new("orbit_type", DataType::Utf8, true)
        ]);

        // Create arrays
        let h_array = Float64Array::from(data.iter().map(|d| d.H).collect::<Vec<_>>());
        let g_array = Float64Array::from(data.iter().map(|d| d.G).collect::<Vec<_>>());
        let epoch_array = Float64Array::from(data.iter().map(|d| d.Epoch).collect::<Vec<_>>());
        let m_array = Float64Array::from(data.iter().map(|d| d.M).collect::<Vec<_>>());
        let peri_array = Float64Array::from(data.iter().map(|d| d.Peri).collect::<Vec<_>>());
        let node_array = Float64Array::from(data.iter().map(|d| d.Node).collect::<Vec<_>>());
        let inc_array = Float64Array::from(data.iter().map(|d| d.i).collect::<Vec<_>>());
        let e_array = Float64Array::from(data.iter().map(|d| d.e).collect::<Vec<_>>());
        let a_array = Float64Array::from(data.iter().map(|d| d.a).collect::<Vec<_>>());
        let principal_desig_array = StringArray::from(
            data.iter().map(|d| d.Principal_desig.as_str()).collect::<Vec<_>>()
        );
        let orbit_type_array = StringArray::from(
            data.iter().map(|d| d.orbit_type.as_deref()).collect::<Vec<_>>()
        );

        // Create RecordBatch
        let batch = RecordBatch::try_new(
            std::sync::Arc::new(schema),
            vec![
                std::sync::Arc::new(h_array),
                std::sync::Arc::new(g_array),
                std::sync::Arc::new(epoch_array),
                std::sync::Arc::new(m_array),
                std::sync::Arc::new(peri_array),
                std::sync::Arc::new(node_array),
                std::sync::Arc::new(inc_array),
                std::sync::Arc::new(e_array),
                std::sync::Arc::new(a_array),
                std::sync::Arc::new(principal_desig_array),
                std::sync::Arc::new(orbit_type_array),
            ],
        ).map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;

        // Write to feather file
        let file = fs::File::create(path)?;
        let mut writer = arrow::ipc::writer::FileWriter::try_new(file, &batch.schema())
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;
        writer.write(&batch)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;
        writer.finish()
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;
        
        Ok(())
    }


    fn read_from_feather(&self, path: &Path) -> PyResult<Vec<MPCData>> {
        let file = fs::File::open(path)?;
        let reader = arrow::ipc::reader::FileReader::try_new(file, None)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;
    
        let mut data = Vec::new();
    
        // Using rayon's parallel iterator for batch processing
        let batches: Vec<_> = reader.collect();
        let processed: Result<Vec<_>, _> = batches.par_iter()
            .map(|batch_result| {
                let batch = batch_result.as_ref()
                    .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;
                self.record_batch_to_mpc_data(batch)
            })
            .collect();
    
        match processed {
            Ok(batch_data) => {
                data.extend(batch_data.into_iter().flatten());
                Ok(data)
            },
            Err(e) => Err(e),
        }
    }

    fn record_batch_to_mpc_data(&self, batch: &RecordBatch) -> PyResult<Vec<MPCData>> {
        let mut data = Vec::with_capacity(batch.num_rows());

        for row in 0..batch.num_rows() {
            let h = batch
                .column_by_name("H")
                .and_then(|col| col.as_any().downcast_ref::<Float64Array>())
                .map(|arr| arr.value(row));

            let g = batch
                .column_by_name("G")
                .and_then(|col| col.as_any().downcast_ref::<Float64Array>())
                .map(|arr| arr.value(row));

            let principal_desig = batch
                .column_by_name("Principal_desig")
                .and_then(|col| col.as_any().downcast_ref::<StringArray>())
                .map(|arr| arr.value(row))
                .unwrap_or_default();

            let orbit_type = batch
                .column_by_name("orbit_type")
                .and_then(|col| col.as_any().downcast_ref::<StringArray>())
                .and_then(|arr| if arr.is_null(row) { None } else { Some(arr.value(row).to_string()) });

            // Get required numeric values
            let get_f64 = |name: &str| -> f64 {
                batch
                    .column_by_name(name)
                    .and_then(|col| col.as_any().downcast_ref::<Float64Array>())
                    .map(|arr| arr.value(row))
                    .unwrap_or_default()
            };

            let mpc_data = MPCData {
                H: h,
                G: g,
                Principal_desig: principal_desig.to_string(),
                orbit_type,
                Epoch: get_f64("Epoch"),
                M: get_f64("M"),
                Peri: get_f64("Peri"),
                Node: get_f64("Node"),
                i: get_f64("i"),
                e: get_f64("e"),
                a: get_f64("a"),
            };

            data.push(mpc_data);
        }

        Ok(data)
    }
}


