use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use std::time::Instant;
use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};

fn main() {
    const N_POINTS: usize = 100000;
    const N_TEST_POINTS: usize = 1000000;
    const BATCH_SIZE: usize = 1000000;
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
        let update = TreeUpdate::new(n_points + i, vertices[i], &delaunay_tree);
        delaunay_tree.insert_point(&update);
    }
    let duration = start.elapsed();
    println!(
        "Time elapsed in constructing the initial tree is: {:?}",
        duration
    );
    println!("Generating test points");
    let mut vertices2: Vec<[f64; 3]> = vec![];
    for _ in 0..N_TEST_POINTS {
        let point = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
        vertices2.push(point);
    }
    println!("Benchmarking update speed");
    // let start = Instant::now();
    // let queue = make_queue(vertices2.clone(), &delaunay_tree);
    // queue
    //     .par_iter()
    //     .enumerate()
    //     .map(|(id, vertex)| TreeUpdate::new(n_points + id, vertex.1, &delaunay_tree))
    //     .collect::<Vec<TreeUpdate<3, 4>>>();
    // let duration = start.elapsed();
    println!("Time elapsed in benchmarking update speed is: {:?}", duration);
    println!("Start inserting test points");
    let start = Instant::now();
    for i in 0..(N_TEST_POINTS / BATCH_SIZE) {
        #[cfg(debug_assertions)]
        {
            let start = Instant::now();
            println!("Starting insert_multiple_points()");
            delaunay_tree
                .add_points_to_tree(vertices2[i * BATCH_SIZE..(i + 1) * BATCH_SIZE].to_vec());
            let duration = start.elapsed();
            println!(
                "Time elapsed in insert_multiple_points() is: {:?}",
                duration
            );
            println!(
                "Number of vertices in the tree: {}",
                delaunay_tree.vertices.len()
            );
        }
        #[cfg(not(debug_assertions))]
        {
            delaunay_tree
                .add_points_to_tree(vertices2[i * BATCH_SIZE..(i + 1) * BATCH_SIZE].to_vec());
        }
    }
    let duration = start.elapsed();
    println!(
        "Time elapsed in inserting test points is: {:?}",
        duration
    );
}
