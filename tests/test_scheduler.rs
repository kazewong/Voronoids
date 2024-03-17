use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};

#[test]
fn test_queuing() {
    let mut vertices = vec![];
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(0.0..1.0);
    for _ in 0..100 {
        let point = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
        vertices.push(point);
    }
    let mut delaunay_tree = DelaunayTree::<3, 4>::new(vertices.clone());
    for i in 0..100 {
        let update = TreeUpdate::new(vertices[i], &delaunay_tree);
        delaunay_tree.insert_point(update);
    }

    let queue = voronoids::scheduler::make_queue(vertices, delaunay_tree);
    println!("{:?}", queue)
}