#![crate_name = "voronoids"]

pub mod delaunay_tree;
pub mod geometry;
pub mod scheduler;

use std::collections::HashMap;

use delaunay_tree::{DelaunayTree, Vertex};
use pyo3::{ffi::PyObject, prelude::*, types::PyDict};

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
}


#[pyfunction]
fn delaunay(points: Vec<[f64; 3]>) -> PyDelauanyTree {
    PyDelauanyTree {
        tree: DelaunayTree::<3, 4>::new(points),
    }
}

#[pymodule]
fn voronoids(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyVertex>()?;
    m.add_class::<PySimplex>()?;
    m.add_class::<PyDelauanyTree>()?;
    m.add_function(wrap_pyfunction!(delaunay, m)?)?;
    Ok(())
}
