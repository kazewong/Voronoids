use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};

#[test]
fn test_delaunay_tree() {
    let vertices = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mut delaunay_tree = DelaunayTree::new(vertices.clone());
    assert_eq!(delaunay_tree.max_simplex_id, 4);
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(0.0..1.0);

    // for i in 0..4 {
    //     println!("{:?}", delaunay_tree.simplices);
    //     println!("{:?}", delaunay_tree.vertices_simplex);
    //     let update = TreeUpdate::new(vertices[i], &delaunay_tree);
    //     delaunay_tree.insert_point(update);
    //     delaunay_tree.check_delaunay();
    // }
    for _ in 0..5 {
        let point = [dist.sample(&mut rng), dist.sample(&mut rng), dist.sample(&mut rng)];
        let update = TreeUpdate::new(point, &delaunay_tree);
        delaunay_tree.insert_point(update);
        delaunay_tree.check_delaunay();
    }
}
