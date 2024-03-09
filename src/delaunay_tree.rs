use std::collections::HashMap;

use kiddo::KdTree;

struct DelaunayTree <T: std::marker::Copy + std::default::Default> {
    kdtree: KdTree<T, 3>,
    vertices: Vec<[T; 3]>,
    vertices_simplex: Vec<Vec<usize>>,

    simplices: HashMap<usize, Vec<usize>>,
    centers: HashMap<usize, [T; 3]>,
    radii: HashMap<usize, T>,
    neighbors: HashMap<usize, Vec<usize>>,
    max_simplex_id: usize,
}

struct TreeUpdate <T>{
    vertex: [T; 3],
    killed_sites: Vec<usize>,
    simplices: Vec<Vec<usize>>,
    simplices_id: Vec<usize>,
    centers: Vec<[T; 3]>,
    radii: Vec<T>,
    neighbors: Vec<(usize, usize)>,
    new_neighbors: Vec<(usize, usize)>,
}

impl <T: std::marker::Copy + std::default::Default> DelaunayTree<T> {
    pub fn new(vertices: Vec<[T; 3]>) -> Self {
        let kdtree = KdTree::new(3);
        let vertices_simplex = vec![vec![]; vertices.len()];
        let vertices = vertices.clone();
        let mut delaunay_tree = DelaunayTree {
            kdtree,
            vertices,
            vertices_simplex,
            simplices: HashMap::new(),
            centers: HashMap::new(),
            radii: HashMap::new(),
            neighbors: HashMap::new(),
            max_simplex_id: 0,
        };
        delaunay_tree.update_tree();
        delaunay_tree
    }

    pub fn locate(&self, vertex: [T; 3]) -> Option<usize> {
        self.kdtree.nearest(&vertex, 1).map(|(id, _)| id)
    }

    pub fn check_delaunay(&self) -> bool{

    }

    pub fn update(&mut self, update: TreeUpdate<T>) {

    }
}

impl <T: std::marker::Copy + std::default::Default> TreeUpdate<T>{
    pub fn new(vertex: [T; 3], tree: DelaunayTree<T>) -> Self {
        TreeUpdate {
            vertex,
            killed_sites: vec![],
            simplices: vec![],
            simplices_id: vec![],
            centers: vec![],
            radii: vec![],
            neighbors: vec![],
            new_neighbors: vec![],
        }
    }
}