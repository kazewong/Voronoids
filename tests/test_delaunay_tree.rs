use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use std::time::Instant;
use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};

#[test]
fn test_delaunay_tree_3d() {
    let mut vertices = vec![];
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(0.0..1.0);
    for _ in 0..1000 {
        let point = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
        vertices.push(point);
    }
    let mut delaunay_tree = DelaunayTree::<3, 4>::new(vertices.clone());
    assert_eq!(delaunay_tree.max_simplex_id, 4);
    let n_points = delaunay_tree.vertices.len();

    for i in 0..100 {
        let update = TreeUpdate::new(n_points+i, vertices[i], &delaunay_tree);
        delaunay_tree.insert_point(&update);
    }
    let mut vertices2: Vec<[f64; 3]> = vec![];
    for _ in 0..1000 {
        let point = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
        vertices2.push(point);
    }
    delaunay_tree.add_points_to_tree(vertices2);
    delaunay_tree.check_delaunay();
}

#[test]
fn test_delaunay_tree_2d() {
    let vertices = vec![[0.3, 0.1], [1.0, 0.2], [0.1, 1.0], [0.5, 0.5]];
    let mut delaunay_tree = DelaunayTree::<2, 3>::new(vertices.clone());
    assert_eq!(delaunay_tree.max_simplex_id, 3);
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(0.0..1.0);
    let n_points = delaunay_tree.vertices.len();


    for i in 0..1000 {
        let start = Instant::now();
        let point = [dist.sample(&mut rng), dist.sample(&mut rng)];
        let update = TreeUpdate::new(n_points+i, point, &delaunay_tree);
        delaunay_tree.insert_point(&update);
        let duration = start.elapsed();
        println!("Point {:?} inserted in {:?}", point, duration);
    }
    delaunay_tree.check_delaunay();
}
