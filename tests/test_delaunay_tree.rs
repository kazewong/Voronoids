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

    // println!("{:?}", delaunay_tree.neighbors);
    // println!("{:?}", delaunay_tree.vertices_simplex);
    // println!("{:?}", delaunay_tree.simplices);

    // let new_point = [0.5, 0.5, 0.5];
    // let update = TreeUpdate::new(new_point, &delaunay_tree);
    // println!("{:?}", update);
    // delaunay_tree.insert_point(update);

    // println!("{:?}", delaunay_tree.neighbors);
    // println!("{:?}", delaunay_tree.vertices_simplex);
    // println!("{:?}", delaunay_tree.simplices);

    for i in 0..100 {
        println!("i: {}", i);
        let new_point = [dist.sample(&mut rng), dist.sample(&mut rng), dist.sample(&mut rng)];
        let update = TreeUpdate::new(new_point, &delaunay_tree);
        delaunay_tree.insert_point(update);
    }
    // assert!(delaunay_tree.check_delaunay());
}