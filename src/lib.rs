#![crate_name = "voronoids"]

pub mod delaunay_tree;
pub mod geometry;
pub mod scheduler;

use std::collections::HashMap;

use delaunay_tree::{DelaunayTree, TreeUpdate};
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
    let mut delaunay_tree = DelaunayTree::<3, 4>::new(points.clone());
    let n_points = delaunay_tree.vertices.len();
    if points.len() > 1e5 as usize {
        println!("More than 1e5 points, using parallel insert");
        for i in 0..1e5 as usize {
            let update = TreeUpdate::new(n_points+i, points[i], &delaunay_tree);
            delaunay_tree.insert_point(&update);
        }
        delaunay_tree.add_points_to_tree(points[1e5 as usize..].to_vec());
    } else {
        println!("Less than 1e5 points, using sequential insert");
        for i in 0..points.len() {
            let update = TreeUpdate::new(n_points+i, points[i], &delaunay_tree);
            delaunay_tree.insert_point(&update);
        }
    }
    PyDelauanyTree {
        tree: delaunay_tree,
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
