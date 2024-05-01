
use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};
use std::time::Instant;

fn main() {
    const N_POINTS: usize = 100000;
    let mut vertices = vec![];
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(0.0..1.0);
    for _ in 0..N_POINTS {
        let point = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
        vertices.push(point);
    }
    let mut delaunay_tree = DelaunayTree::<3, 4>::new(vertices.clone());
    let n_points = delaunay_tree.vertices.len();
    println!("Start constructing the initial tree");
    let start = Instant::now();
    for i in 0..N_POINTS {
        let update = TreeUpdate::new(n_points+i, vertices[i], &delaunay_tree);
        delaunay_tree.insert_point(&update);
    }
    let duration = start.elapsed();
    println!("Number of vertices in the initial tree: {}", delaunay_tree.vertices.len());
    println!("Time elapsed in constructing the initial tree is: {:?}", duration);
}