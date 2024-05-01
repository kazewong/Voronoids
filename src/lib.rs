#![crate_name = "voronoids"]

pub mod delaunay_tree;
pub mod geometry;
pub mod scheduler;

use delaunay_tree::DelaunayTree;
use pyo3::prelude::*;

#[pyclass]
struct PyDelauanyTree {
    tree: DelaunayTree<3, 4>,
}

#[pyfunction]
fn delaunay(points: Vec<[f64; 3]>)-> PyDelauanyTree {
    PyDelauanyTree {
        tree: DelaunayTree::<3, 4>::new(points),
    }
}
