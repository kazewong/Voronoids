use crate::geometry::{circumsphere, in_sphere};
use kiddo::{KdTree, SquaredEuclidean};
use nalgebra::Point3;
use parry3d_f64::bounding_volume::details::point_cloud_bounding_sphere;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct DelaunayTree {
    pub kdtree: KdTree<f64, 3>,
    pub vertices: Vec<[f64; 3]>,
    pub vertices_simplex: Vec<Vec<usize>>,
    pub simplices: HashMap<usize, [usize; 4]>,
    pub centers: HashMap<usize, [f64; 3]>,
    pub radii: HashMap<usize, f64>,
    pub neighbors: HashMap<usize, Vec<usize>>,
    pub max_simplex_id: usize,
}

pub struct TreeUpdate {
    vertex: [f64; 3],
    killed_sites: Vec<usize>,
    simplices: Vec<[usize;4]>,
    simplices_id: Vec<usize>,
    centers: Vec<[f64; 3]>,
    radii: Vec<f64>,
    neighbors: Vec<(usize, usize)>,
    new_neighbors: Vec<(usize, usize)>,
}

impl DelaunayTree {
    pub fn new(vertices: Vec<[f64; 3]>) -> Self {
        // Turn vertices into nalgebra points
        let points: Vec<Point3<f64>> = vertices
            .iter()
            .map(|v| Point3::new(v[0], v[1], v[2]))
            .collect();
        let (center, mut radius) = point_cloud_bounding_sphere(&points);
        radius *= 5.0;

        let first_vertex = [0. + center[0], 0. + center[1], radius + center[2]];
        let second_vertex = [radius + center[0], 0. + center[1], -radius + center[2]];
        let third_vertex = [
            -radius / 2. + center[0],
            radius * 3f64.sqrt() / 2. + center[1],
            -radius + center[2],
        ];
        let fourth_vertex = [
            -radius / 2. + center[0],
            -radius * 3f64.sqrt() / 2. + center[1],
            -radius + center[2],
        ];
        let ghost_vertex = [
            [0. + center[0], 0. + center[1], radius + center[2]],
            [0. + center[0], 0. + center[1], radius + center[2]],
            [0. + center[0], 0. + center[1], radius + center[2]],
            [radius + center[0], 0. + center[1], -radius + center[2]],
        ];
        let mut kdtree = KdTree::new();

        kdtree.add(&first_vertex, 0);
        kdtree.add(&second_vertex, 1);
        kdtree.add(&third_vertex, 2);
        kdtree.add(&fourth_vertex, 3);
        for (i, vertex) in ghost_vertex.iter().enumerate() {
            kdtree.add(vertex, (i + 4) as u64);
        }

        let vertices_simplex = [
            vec![0, 1, 2, 3],
            vec![0, 1, 3, 4],
            vec![0, 1, 2, 4],
            vec![0, 2, 3, 4],
            vec![1],
            vec![2],
            vec![3],
            vec![4],
        ]
        .to_vec();
        let simplices = HashMap::from([
            (0, [0, 1, 2, 3]),
            (1, [4, 0, 1, 2]),
            (2, [5, 0, 2, 3]),
            (3, [6, 0, 3, 1]),
            (4, [7, 1, 2, 3]),
        ]);
        let centers = HashMap::from([
            (0, center.into()),
            (1, [0., 0., 0.]),
            (2, [0., 0., 0.]),
            (3, [0., 0., 0.]),
            (4, [0., 0., 0.]),
        ]);
        let radii = HashMap::from([(0, radius), (1, 0.), (2, 0.), (3, 0.), (4, 0.)]);
        let neighbors = HashMap::from([
            (0, vec![1, 2, 3, 4]),
            (1, vec![0]),
            (2, vec![0]),
            (3, vec![0]),
            (4, vec![0]),
        ]);
        let delaunay_tree = DelaunayTree {
            kdtree,
            vertices,
            vertices_simplex,
            simplices,
            centers,
            radii,
            neighbors,
            max_simplex_id: 4,
        };
        delaunay_tree
    }

    pub fn check_delaunay(&self) -> bool {
        for (id, simplex) in self.simplices.iter() {
            for vertex in self.vertices.iter() {
                if in_sphere(
                    *vertex,
                    *self.centers.get(id).unwrap(),
                    *self.radii.get(id).unwrap(),
                ) && !simplex.contains(&id)
                {
                    return false;
                }
            }
        }
        true
    }

    pub fn insert_point(&mut self, update: TreeUpdate) {
        let mut killed_sites = update.killed_sites;
        killed_sites.sort();
        self.vertices.push(update.vertex);

        for (i, (neighbor_id, killed_id)) in update.neighbors.iter().enumerate() {
            for j in 0..self.neighbors[neighbor_id].len() {
                if self.neighbors[neighbor_id][j] == *killed_id {
                    self.neighbors.get_mut(neighbor_id).unwrap()[j] =
                        self.max_simplex_id + update.simplices_id[i];
                }
            }
        }
    }

    pub fn get_new_simplices(
        &self,
        killed_site_id: usize,
        vertex: [f64; 3],
        vertex_id: usize,
    ) -> (
        Vec<Vec<usize>>,
        Vec<[f64; 3]>,
        Vec<f64>,
        Vec<(usize, usize)>,
    ) {
        let mut simplices: Vec<Vec<usize>> = vec![];
        let mut simplices_id: Vec<usize> = vec![];
        let mut centers: Vec<[f64; 3]> = vec![];
        let mut radii: Vec<f64> = vec![];
        let mut neighbors: Vec<(usize, usize)> = vec![];

        let killed_site: [usize; 4] = *self.simplices.get(&killed_site_id).unwrap();

        for neighbor_id in self.neighbors.get(&killed_site_id).unwrap() {
            let neighbor_simplex = *self.simplices.get(neighbor_id).unwrap();
            if !in_sphere(vertex, *self.centers.get(neighbor_id).unwrap(), *self.radii.get(neighbor_id).unwrap()) {
                let mut new_simplex = [0, 0, 0, 0];
                for i in 0..4 {
                    if killed_site.contains(&neighbor_simplex[i]) {
                        new_simplex[i] = neighbor_simplex[i];
                    } else {
                        new_simplex[i] = vertex_id;
                    }
                }
                let (center, radius) = circumsphere([
                    self.vertices[new_simplex[0]],
                    self.vertices[new_simplex[1]],
                    self.vertices[new_simplex[2]],
                    self.vertices[new_simplex[3]],
                ]);
                simplices.push(new_simplex.to_vec());
                simplices_id.push(self.max_simplex_id + simplices.len());
                centers.push(center);
                radii.push(radius);
                neighbors.push((*neighbor_id, killed_site_id));
            }
        }
        (simplices, centers, radii, neighbors)
    }

    fn find_all_neighbors(
        self,
        output: &mut Vec<usize>,
        node_id: usize,
        vertex: [f64; 3],
    ) -> Vec<usize> {
        let neighbors = self.neighbors.get(&node_id).unwrap();
        for neighbor in neighbors {
            if !output.contains(neighbor)
                && in_sphere(
                    vertex,
                    *self.centers.get(neighbor).unwrap(),
                    *self.radii.get(neighbor).unwrap(),
                )
            {
                output.push(*neighbor);
                self.clone().find_all_neighbors(output, *neighbor, vertex);
            }
        }
        output.to_vec()
    }

    fn locate(&self, vertex: [f64; 3]) -> Vec<usize> {
        let mut output: Vec<usize> = vec![];
        let simplex_id = &self.vertices_simplex
            [self.kdtree.nearest_one::<SquaredEuclidean>(&vertex).item as usize];
        for id in simplex_id {
            if in_sphere(
                vertex,
                *self.centers.get(id).unwrap(),
                *self.radii.get(id).unwrap(),
            ) {
                output.push(*id);
                output = self.clone().find_all_neighbors(&mut output, *id, vertex);
            }
        }
        output
    }
}

impl TreeUpdate {
    pub fn new(vertex: [f64; 3], tree: DelaunayTree) -> Self {
        let killed_sites = tree.locate(vertex);
        let id = tree.vertices.len() + 1;
        let mut simplices: Vec<Vec<usize>> = vec![];
        let mut simplices_id: Vec<usize> = vec![];
        let mut centers: Vec<[f64; 3]> = vec![];
        let mut radii: Vec<f64> = vec![];
        let mut neighbors: Vec<(usize, usize)> = vec![];

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
