use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};
use rand::prelude::*;
use rand::distributions::{Uniform, Distribution};

#[test]
fn test_delaunay_tree() {
    let vertices = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mut delaunay_tree = DelaunayTree::new(vertices);
    assert_eq!(delaunay_tree.max_simplex_id, 4);
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(0.0..1.0);

    for i in 0..10 {
        println!("i: {}", i);
        let new_point = [dist.sample(&mut rng), dist.sample(&mut rng), dist.sample(&mut rng)];
        let update = TreeUpdate::new(new_point, &delaunay_tree);
        delaunay_tree.insert_point(update);
    // println!("{:?}", delaunay_tree.neighbors);
    // println!("{:?}", delaunay_tree.vertices_simplex);
    }
    assert!(delaunay_tree.check_delaunay());
}