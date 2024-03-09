use std::collections::HashMap;
use kiddo::KdTree;
use nalgebra::Point3;
use parry3d_f64::bounding_volume::details::point_cloud_bounding_sphere;
use crate::geometry::circumsphere;

#[derive(Debug, Clone)]
struct DelaunayTree{
    kdtree: KdTree<f64, 3>,
    vertices: Vec<[f64; 3]>,
    vertices_simplex: Vec<Vec<usize>>,

    simplices: HashMap<usize, Vec<usize>>,
    centers: HashMap<usize, [f64; 3]>,
    radii: HashMap<usize, f64>,
    neighbors: HashMap<usize, Vec<usize>>,
    max_simplex_id: usize,
}

struct TreeUpdate{
    vertex: [f64; 3],
    killed_sites: Vec<usize>,
    simplices: Vec<Vec<usize>>,
    simplices_id: Vec<usize>,
    centers: Vec<[f64; 3]>,
    radii: Vec<f64>,
    neighbors: Vec<(usize, usize)>,
    new_neighbors: Vec<(usize, usize)>,
}

impl DelaunayTree{
    pub fn new(vertices: Vec<[f64; 3]>) -> Self {
        // Turn vertices into nalgebra points
        let points: Vec<Point3<f64>> = vertices.iter().map(|v| Point3::new(v[0], v[1], v[2])).collect();
        let (center, mut radius) = point_cloud_bounding_sphere(&points);
        radius *= 5.0;


        let first_vertex = [0.+center[0], 0. + center[1], radius+center[2]];
        let second_vertex = [radius+center[0], 0. + center[1], -radius+center[2]];
        let third_vertex = [-radius / 2.+center[0], radius * 3f64.sqrt() / 2.+center[1], -radius+center[2]];
        let fourth_vertex = [-radius / 2.+center[0], -radius * 3f64.sqrt() / 2. +center[1], -radius+center[2]];
        let ghost_vertex = [[0.+center[0], 0. + center[1], radius+center[2]], [0.+center[0], 0. + center[1], radius+center[2]], [0.+center[0], 0. + center[1], radius+center[2]], [radius+center[0], 0. + center[1], -radius+center[2]]];
        let mut kdtree = KdTree::new();

        kdtree.add(&first_vertex, 0);
        kdtree.add(&second_vertex, 1);
        kdtree.add(&third_vertex, 2);
        kdtree.add(&fourth_vertex, 3);
        for (i, vertex) in ghost_vertex.iter().enumerate() {
            kdtree.add(vertex, (i + 4) as u64);
        }

        let vertices_simplex = [vec![1, 2, 3, 4], vec![1, 2, 4, 5], vec![1, 2, 3, 5], vec![1, 3, 4, 5], vec![2], vec![3], vec![4], vec![5]].to_vec();
        let simplices = HashMap::from([(1, vec![1, 2, 3, 4]), (2, vec![5, 1, 2, 3]), (3, vec![6, 1, 3, 4]), (4, vec![7, 1, 4, 2]), (5, vec![8, 2, 3, 4])]);
        let centers = HashMap::from([(1, center.into()), (2, [0., 0., 0.]), (3, [0., 0., 0.]), (4, [0., 0., 0.]), (5, [0., 0., 0.])]);
        let radii = HashMap::from([(1, radius), (2, 0.), (3, 0.), (4, 0.), (5, 0.)]);
        let neighbors = HashMap::from([(1, vec![2, 3, 4, 5]), (2, vec![1]), (3, vec![1]), (4, vec![1]), (5, vec![1])]);
        let delaunay_tree = DelaunayTree {
            kdtree,
            vertices,
            vertices_simplex,
            simplices,
            centers,
            radii,
            neighbors,
            max_simplex_id: 5,
        };
        delaunay_tree
    }

    // pub fn locate(&self, vertex: [T; 3]) -> Option<usize> {
    //     self.kdtree.nearest(&vertex, 1).map(|(id, _)| id)
    // }

    // pub fn check_delaunay(&self) -> bool{

    // }

    // pub fn update(&mut self, update: TreeUpdate<T>) {

    // }
}

// impl <T: std::marker::Copy + std::default::Default> TreeUpdate<T>{
//     pub fn new(vertex: [T; 3], tree: DelaunayTree<T>) -> Self {
//         TreeUpdate {
//             vertex,
//             killed_sites: vec![],
//             simplices: vec![],
//             simplices_id: vec![],
//             centers: vec![],
//             radii: vec![],
//             neighbors: vec![],
//             new_neighbors: vec![],
//         }
//     }
// }