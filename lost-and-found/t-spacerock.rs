// use pyo3::prelude::*;
// use pyo3::types::PyType;

// use reqwest::blocking::Client;
// use std::collections::HashMap;
// use serde_json;

// #[pyclass]
// #[derive(Clone)]
// pub struct SpaceRock {
//     #[pyo3(get, set)]
//     pub name: String,
//     #[pyo3(get, set)]
//     pub position: Vec<f64>,
//     #[pyo3(get, set)]
//     pub velocity: Vec<f64>,
// }

// #[pymethods]
// impl SpaceRock {
//     #[new]
//     fn new(name: String) -> Self {
//         SpaceRock { name, position: vec![0.0, 0.0, 0.0], velocity: vec![0.0, 0.0, 0.0] }
//     }

//     fn __repr__(&self) -> String {
//         format!("SpaceRock: {}", self.name)
//     }

//     #[classmethod]
//     fn from_state(cls: &PyType, name: &str, position: (f64, f64, f64), velocity: (f64, f64, f64)) -> Self {
//         SpaceRock { name: name.to_string(), position: vec![position.0, position.1, position.2], velocity: vec![velocity.0, velocity.1, velocity.2] }
//     }

//     #[classmethod]
//     fn from_horizons(cls: &PyType, name: &str) -> Self {
//         let client = reqwest::blocking::Client::new();

//         let mut params = HashMap::new();
//         // make the command the name, but with a single quote around it as a string

//         let command_str = format!("'{}'", name);

//         params.insert("command", command_str.as_str());
//         params.insert("start_time", "'2023-03-14'");
//         params.insert("stop_time", "'2023-03-15'");
//         params.insert("make_ephem", "'yes'");
//         params.insert("ephem_type", "'vectors'");
//         params.insert("center", "'@ssb'");
//         params.insert("ref_plane", "'ecliptic'");
//         params.insert("step_size", "'1d'");
//         params.insert("ref_system", "'J2000'");
//         params.insert("vec_corr", "'None'");
//         params.insert("out_units", "'AU-D'");
//         params.insert("csv_format", "'yes'");
//         params.insert("vec_delta_t", "'no'");
//         params.insert("vec_table", "'2x'");
//         params.insert("vec_labels", "'no'");

    
//         let response = client.get("https://ssd.jpl.nasa.gov/api/horizons.api?")
//             .query(&params)
//             .send().unwrap();

//         let json: serde_json::Value = response.json().unwrap();
//         let text = json["result"].as_str();

//         let lines: Vec<&str> = text.ok_or_else(| | "No data").unwrap().split('\n').collect();
//         let first_data_line = lines.iter().skip_while(|&line| !line.starts_with("$$SOE")).skip(1).next().ok_or("No data").unwrap();
        
//         let data: Vec<f64> = first_data_line.split(',').filter_map(|s| s.trim().parse::<f64>().ok()).collect();
//         // let epoch = Time { epoch: data[0], timescale: "tdb".to_string(), format: "jd".to_string() };
//         let (x, y, z, vx, vy, vz) = (data[1], data[2], data[3], data[4], data[5], data[6]);

//         SpaceRock { name: name.to_string(), position: vec![x, y, z], velocity: vec![vx, vy, vz] }

//     }

//     pub fn hvec(&self) -> (f64, f64, f64) {
//         // specific angular momentum vector
//         let hx = self.position[1]*self.velocity[2] - self.position[2]*self.velocity[1];
//         let hy = self.position[2]*self.velocity[0] - self.position[0]*self.velocity[2];
//         let hz = self.position[0]*self.velocity[1] - self.position[1]*self.velocity[0];
//         (hx, hy, hz)
//     }

//     pub fn h(&self) -> f64 {
//         let (hx, hy, hz) = self.hvec();
//         (hx*hx + hy*hy + hz*hz).sqrt()
//     }


//     pub fn rsq(&self) -> f64 {
//         // self.position.iter().map(|x| x*x).sum()
//         self.position[0]*self.position[0] + self.position[1]*self.position[1] + self.position[2]*self.position[2]
//     }

//     pub fn vsq(&self) -> f64 {
//         self.velocity.iter().map(|x| x*x).sum()
//     }

//     // implement a clone method
//     pub fn copy(&self) -> Self {
//         SpaceRock { name: self.name.clone(), position: self.position.clone(), velocity: self.velocity.clone() }
//     }
// }
