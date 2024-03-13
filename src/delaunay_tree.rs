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

#[derive(Debug, Clone)]
pub struct TreeUpdate {
    vertex: [f64; 3],
    killed_sites: Vec<usize>,
    simplices: Vec<[usize; 4]>,
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

        let vertices = vec![
            first_vertex,
            second_vertex,
            third_vertex,
            fourth_vertex,
            ghost_vertex[0],
            ghost_vertex[1],
            ghost_vertex[2],
            ghost_vertex[3],
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
                ) && !simplex.contains(&id) && simplex.iter().all(|&x| {
                    x > 7
                })
                {
                    return false;
                }
            }
        }
        true
    }

    pub fn insert_point(&mut self, update: TreeUpdate) {
        let mut killed_sites = update.killed_sites;
        self.kdtree.add(&update.vertex, self.vertices.len() as u64);
        killed_sites.sort();
        self.vertices.push(update.vertex);

        // Update simplices
        for (i, simplex) in update.simplices.iter().enumerate() {
            let current_id = self.max_simplex_id + update.simplices_id[i];
            self.simplices.insert(current_id, *simplex);
            self.centers.insert(current_id, update.centers[i]);
            self.radii.insert(current_id, update.radii[i]);
            self.neighbors.insert(
                current_id,
                vec![update.neighbors[i].0],
            );
        }

        // update neighbor relations

        for (i, (neighbor_id, killed_id)) in update.neighbors.iter().enumerate() {
            for j in 0..self.neighbors[neighbor_id].len() {
                if self.neighbors[neighbor_id][j] == *killed_id {
                    self.neighbors.get_mut(neighbor_id).unwrap()[j] =
                        self.max_simplex_id + update.simplices_id[i];
                }
            }
        }

        for (new_neighbor_id1, new_neighbor_id2) in update.new_neighbors.iter(){
            self.neighbors.get_mut(&(self.max_simplex_id+ *new_neighbor_id1)).unwrap().push(self.max_simplex_id + *new_neighbor_id2);
        }

        // Update vertices_simplex

        self.vertices_simplex.push(vec![]);
        for i in 0..update.simplices.len() {
            for j in 0..4 {
                self.vertices_simplex[update.simplices[i][j]].push(self.max_simplex_id + update.simplices_id[i]);
            }
        }
        for killed_site_id in killed_sites.iter() {
            for i in 0..4 {
                self.vertices_simplex[self.simplices[killed_site_id][i]].retain(|&x| x != *killed_site_id);
            }
        }

        // Remove killed sites
        for killed_sites_id in killed_sites.iter() {
            self.simplices.remove(killed_sites_id);
            self.centers.remove(killed_sites_id);
            self.radii.remove(killed_sites_id);
            self.neighbors.remove(killed_sites_id);
        }

        self.max_simplex_id += update.simplices.len();
    }

    pub fn get_new_simplices(
        &self,
        killed_site_id: usize,
        vertex: [f64; 3],
        vertex_id: usize,
    ) -> (
        Vec<[usize; 4]>,
        Vec<[f64; 3]>,
        Vec<f64>,
        Vec<(usize, usize)>,
    ) {
        let mut simplices: Vec<[usize; 4]> = vec![];
        let mut simplices_id: Vec<usize> = vec![];
        let mut centers: Vec<[f64; 3]> = vec![];
        let mut radii: Vec<f64> = vec![];
        let mut neighbors: Vec<(usize, usize)> = vec![];

        let killed_site: [usize; 4] = *self.simplices.get(&killed_site_id).unwrap();

        for neighbor_id in self.neighbors.get(&killed_site_id).unwrap() {
            let neighbor_simplex = *self.simplices.get(neighbor_id).unwrap();
            if !in_sphere(
                vertex,
                *self.centers.get(neighbor_id).unwrap(),
                *self.radii.get(neighbor_id).unwrap(),
            ) {
                let mut new_simplex = [vertex_id, 0, 0, 0];
                let mut count = 1;
                for i in 0..4 {
                    if killed_site.contains(&neighbor_simplex[i]) {
                        new_simplex[count] = neighbor_simplex[i];
                        count += 1;
                    }
                }
                let (center, radius) = circumsphere([
                    vertex,
                    self.vertices[new_simplex[1]],
                    self.vertices[new_simplex[2]],
                    self.vertices[new_simplex[3]],
                ]);
                simplices.push(new_simplex);
                simplices_id.push(self.max_simplex_id + simplices.len());
                centers.push(center);
                radii.push(radius);
                neighbors.push((*neighbor_id, killed_site_id));
            }
        }
        (simplices, centers, radii, neighbors)
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
                if !output.contains(id) {
                    output.push(*id);
                }                    
                output = self.clone().find_all_neighbors(&mut output, *id, vertex);
            }
        }
        output
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
}

fn pair_simplices(simplices: &Vec<[usize; 4]>, simplices_id: &Vec<usize>) -> Vec<(usize, usize)> {
    let mut new_neighbors: Vec<(usize, usize)> = vec![];
    let n_simplices = simplices.len();
    for i in 0..n_simplices {
        for j in i..n_simplices {
            if i != j {
                let mut count = 0;
                for k in 0..4 {
                    if simplices[i].contains(&simplices[j][k]) {
                        count += 1;
                    }
                }
                if count == 3 {
                    new_neighbors.push((simplices_id[i], simplices_id[j]));
                    new_neighbors.push((simplices_id[j], simplices_id[i]));
                }
            }
        }
    }
    new_neighbors
}

impl TreeUpdate {
    pub fn new(vertex: [f64; 3], tree: &DelaunayTree) -> Self {
        let killed_sites = tree.locate(vertex);
        let id = tree.vertices.len();
        let mut simplices: Vec<[usize; 4]> = vec![];
        let mut simplices_id: Vec<usize> = vec![];
        let mut centers: Vec<[f64; 3]> = vec![];
        let mut radii: Vec<f64> = vec![];
        let mut neighbors: Vec<(usize, usize)> = vec![];

        let mut simplices_counter = 1;

        for i in 0..killed_sites.len() {
            let (simplices_, centers_, radii_, neighbors_) =
                tree.get_new_simplices(killed_sites[i], vertex, id);
            simplices.extend(simplices_.clone());
            centers.extend(centers_);
            radii.extend(radii_);
            neighbors.extend(neighbors_);
            simplices_id.extend(
                (simplices_counter..simplices_counter + simplices_.len()).collect::<Vec<usize>>(),
            );
            simplices_counter += simplices_.len();
        }

        let new_neighbors: Vec<(usize, usize)> = pair_simplices(&simplices, &simplices_id);

        TreeUpdate {
            vertex,
            killed_sites,
            simplices,
            simplices_id,
            centers,
            radii,
            neighbors,
            new_neighbors,
        }
    }
}
