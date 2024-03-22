use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};
use voronoids::scheduler::find_placement;

#[test]
fn test_queuing() {
    const N_POINTS: usize = 1000;
    let mut vertices = vec![];
    let mut vertices2 = vec![];
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(0.0..1.0);
    for _ in 0..N_POINTS {
        let point = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
        vertices.push(point);
        let point2 = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
        vertices2.push(point2);
    }
    let mut delaunay_tree = DelaunayTree::<3, 4>::new(vertices.clone());
    for i in 0..N_POINTS {
        let update = TreeUpdate::new(vertices[i], &delaunay_tree);
        delaunay_tree.insert_point(update);
    }

    let queue = voronoids::scheduler::make_queue(vertices2, &delaunay_tree);
<<<<<<< HEAD
    let placement = find_placement(&queue);
=======
    let placement = find_placement(queue);
>>>>>>> d98f7fc (make_queue should use reference)
    println!("{:?}", placement);
}