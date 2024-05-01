#![crate_name = "voronoids"]

pub mod delaunay_tree;
pub mod geometry;
pub mod scheduler;

use std::collections::HashMap;

use delaunay_tree::{DelaunayTree};
use pyo3::prelude::*;

#[pyclass]
struct PyVertex {
    point: [f64; 3],
    simplex: Vec<usize>,
}

#[pymethods]
impl PyVertex {
    #[getter]
    fn point(&self) -> [f64; 3] {
        self.point
    }

    #[getter]
    fn simplex(&self) -> Vec<usize> {
        self.simplex.clone()
    }
}

#[pyclass]
struct PySimplex {
    vertices: [usize; 4],
    center: [f64; 3],
    radius: f64,
    neighbors: Vec<usize>,
}

#[pymethods]
impl PySimplex {
    #[getter]
    fn vertices(&self) -> [usize; 4] {
        self.vertices
    }

    #[getter]
    fn center(&self) -> [f64; 3] {
        self.center
    }

    #[getter]
    fn radius(&self) -> f64 {
        self.radius
    }

    #[getter]
    fn neighbors(&self) -> Vec<usize> {
        self.neighbors.clone()
    }
}

#[pyclass]
struct PyDelauanyTree {
    tree: DelaunayTree<3, 4>,
}
#[pymethods]
impl PyDelauanyTree {
    #[getter]
    fn max_simplex_id(&self) -> usize {
        self.tree.max_simplex_id
    }

    #[getter]
    fn vertices(&self) -> HashMap<usize, PyVertex> {
        HashMap::from_iter(self.tree.vertices.iter().map(|id| {
            let vertex = &self.tree.vertices.get(id.key()).unwrap();
            (
                id.key().clone(),
                PyVertex {
                    point: vertex.coordinates,
                    simplex: vertex.simplex.clone(),
                },
            )
        }))
    }

    #[getter]
    fn simplices(&self) -> HashMap<usize, PySimplex> {
        HashMap::from_iter(self.tree.simplices.iter().map(|id| {
            let simplex = &self.tree.simplices.get(id.key()).unwrap();
            (
                id.key().clone(),
                PySimplex {
                    vertices: simplex.vertices,
                    center: simplex.center,
                    radius: simplex.radius,
                    neighbors: simplex.neighbors.clone(),
                },
            )
        }))
    }
}

#[pyfunction]
fn delaunay(points: Vec<[f64; 3]>) -> PyDelauanyTree {
    PyDelauanyTree {
        tree: DelaunayTree::<3, 4>::new(points),
    }
}

#[pymodule]
fn voronoids(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyVertex>()?;
    m.add_class::<PySimplex>()?;
    m.add_class::<PyDelauanyTree>()?;
    m.add_function(wrap_pyfunction!(delaunay, m)?)?;
    Ok(())
}
