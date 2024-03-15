use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};

#[test]
fn test_delaunay_tree3D() {
    let vertices = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mut delaunay_tree = DelaunayTree::<3,4>::new(vertices.clone());
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
    for i in 0..10 {
        println!("{:?}", i);
        let point = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
        let update = TreeUpdate::new(point, &delaunay_tree);
        let mut keys = delaunay_tree.simplices.keys().collect::<Vec<_>>();
        keys.sort();
        for key in keys {
            println!(
                "{:?}, {:?}, {:?}",
                &key, delaunay_tree.simplices[&key], delaunay_tree.neighbors[&key]
            );
        }
        delaunay_tree.insert_point(update);
        delaunay_tree.check_delaunay();
    }
}

#[test]
fn test_delaunay_tree2D() {
    let vertices = vec![
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
    ];
    let mut delaunay_tree = DelaunayTree::<2,3>::new(vertices.clone());
    assert_eq!(delaunay_tree.max_simplex_id, 3);
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(0.0..1.0);

    for i in 0..4 {
        println!("{:?}", delaunay_tree.simplices);
        println!("{:?}", delaunay_tree.vertices);
        let update = TreeUpdate::new(vertices[i], &delaunay_tree);
        println!("{:?}", update);
        delaunay_tree.insert_point(update);
        delaunay_tree.check_delaunay();
    }
    // for i in 0..10 {
    //     println!("{:?}", i);
    //     let point = [
    //         dist.sample(&mut rng),
    //         dist.sample(&mut rng),
    //         dist.sample(&mut rng),
    //     ];
    //     let update = TreeUpdate::new(point, &delaunay_tree);
    //     let mut keys = delaunay_tree.simplices.keys().collect::<Vec<_>>();
    //     keys.sort();
    //     for key in keys {
    //         println!(
    //             "{:?}, {:?}, {:?}",
    //             &key, delaunay_tree.simplices[&key], delaunay_tree.neighbors[&key]
    //         );
    //     }
    //     delaunay_tree.insert_point(update);
    //     delaunay_tree.check_delaunay();
    // }
}
