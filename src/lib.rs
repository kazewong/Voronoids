#![crate_name = "voronoids"]

pub mod delaunay_tree;
pub mod geometry;
pub mod scheduler;

use dashmap::DashMap;
use delaunay_tree::{DelaunayTree, Vertex};
use pyo3::{ffi::PyObject, prelude::*, types::PyDict};

#[pyclass]
struct PyDelauanyTree {
    tree: DelaunayTree<3, 4>,
}

#[pyclass]
struct PyVertex {
    point: [f64; 3],
    neighbors: Vec<usize>,
}

#[pyclass]
struct PySimplex {
    vertices: [usize; 4],
    center: [f64; 3],
    radius: f64,
    neighbors: Vec<usize>,
}

#[pymethods]
impl PyDelauanyTree{
    #[getter]
    fn max_simplex_id(&self) -> usize {
        self.tree.max_simplex_id
    }

}

#[pyfunction]
fn delaunay(points: Vec<[f64; 3]>)-> PyDelauanyTree {
    PyDelauanyTree {
        tree: DelaunayTree::<3, 4>::new(points),
    }
}

#[pymodule]
fn voronoids(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyDelauanyTree>()?;
    m.add_function(wrap_pyfunction!(delaunay, m)?)?;
    Ok(())
}
