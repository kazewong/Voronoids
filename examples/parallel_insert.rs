use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};
use std::time::Instant;

fn main() {
    const N_POINTS: usize = 10000;
    const N_TEST_POINTS: usize = 10000;
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
        delaunay_tree.insert_point(update);
    }
    let duration = start.elapsed();
    println!("Time elapsed in constructing the initial tree is: {:?}", duration);

    for i in 0..10{
        let mut vertices2: Vec<[f64; 3]> = vec![];
        for _ in 0..N_TEST_POINTS {
            let point = [
                dist.sample(&mut rng),
                dist.sample(&mut rng),
                dist.sample(&mut rng),
            ];
            vertices2.push(point);
        }
        let start = Instant::now();
        println!("Starting insert_multiple_points()");
        delaunay_tree.insert_multiple_points(vertices2);
        let duration = start.elapsed();
        println!("Time elapsed in insert_multiple_points() is: {:?}", duration);
        println!("Number of vertices in the tree: {}", delaunay_tree.vertices.len());
    }
}