// // pyclass that holds references to a list of SpaceRock objects
// use pyo3::prelude::*;
// use pyo3::exceptions::PyIndexError;
// use pyo3::types::PyDict;
// use pyo3::types::PyTuple;

// use crate::spacerock::spacerock::SpaceRock;

// use tabled::{
//     builder::Builder,
//     settings::{style::Style, themes::Colorization, Color},
// };
// use tabled::settings::Padding;
// use tabled::Table;
// use tabled::settings::Shadow;



// #[pyclass]
// pub struct RockCollection {
//     #[pyo3(get, set)]
//     pub rocks: Vec<Py<SpaceRock>>,
// }

// #[pymethods]
// impl RockCollection {
//     #[new]
//     fn new() -> Self {
//         RockCollection { rocks: Vec::new() }
//     }

//     #[staticmethod]
//     fn from_rocks(rocks: Vec<Py<SpaceRock>>) -> Self {
//         RockCollection { rocks }
//     }

//     fn add(&mut self, rock: Py<SpaceRock>) {
//         self.rocks.push(rock);
//     }


//     fn h(&self, py: Python) -> Vec<f64> {
//         // return a vector of h values
//         let mut h = Vec::new();
//         for rock in self.rocks.iter() {
//             let rock = rock.borrow(py);
//             h.push(rock.h());
//         }
//         h
//     }

//     // make indexable
//     fn __getitem__(&self, index: usize) -> PyResult<Py<SpaceRock>> {
//         if index < self.rocks.len() {
//             Ok(self.rocks[index].clone())
//         } else {
//             Err(PyIndexError::new_err("Index out of range!"))
//         }
//     }

//     fn __len__(&self) -> usize {
//         self.rocks.len()
//     }

//     fn keep(&mut self, py: Python, indices: Vec<bool>) -> Self {
//         if indices.len() != self.rocks.len() {
//             panic!("Length of indices must match length of rocks!");
//         }

//         let mut new_rocks = Vec::new();
//         for (i, rock) in self.rocks.iter().enumerate() {
//             if indices[i] {
//                 new_rocks.push(rock.clone());
//             }
//         }
//         RockCollection { rocks: new_rocks }
//     }

//     fn copy(&self) -> Self {
//         RockCollection { rocks: self.rocks.clone() }
//     }

//     fn __repr__(&self, py: Python) -> String {

//         let mut data = Vec::new();
//         data.push(vec!["Name".to_string(), "x (au)".to_string(), "y (au)".to_string(), "z (au)".to_string()]);
//         for rock in self.rocks.iter() {
//             let rock = rock.borrow(py);
//             // data.push(vec![&rock.name.clone(), &rock.position[0].clone().to_string()]);
//             let name = rock.name.clone();
//             let x = format!("{:.3}", rock.position.clone()[0]);
//             let y = format!("{:.3}", rock.position.clone()[1]);
//             let z = format!("{:.3}", rock.position.clone()[2]);
//             let row = vec![name, x, y, z];
//             data.push(row);
//         }
//         let mut table = Table::from_iter(data);
//         // table.with(Style::rounded());
//         // table.with(Padding::new(2, 2, 0, 0));
//         table.with(Style::modern());
//         // table.with(Shadow::new(3).set_offset(6));
//         table.to_string()
//     }


//     // fn perturb_rocks(&mut self, py: Python) {
//     //     for rock in self.rocks.iter() {
//     //         let mut borrow = rock.borrow_mut(py);
//     //         borrow.position[0] += 0.1;
//     //         borrow.position[1] += 0.1;
//     //         borrow.position[2] += 0.1;
//     //         borrow.velocity[0] += 0.1;
//     //         borrow.velocity[1] += 0.1;
//     //         borrow.velocity[2] += 0.1;
//     //     }
//     // }
// }
