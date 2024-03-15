use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};
use plotters::prelude::*;

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

const OUT_FILE_NAME: &str = "./test.png";

#[test]
fn test_delaunay_tree2D() {
    let root = BitMapBackend::new(OUT_FILE_NAME, (1024, 768)).into_drawing_area();

    root.fill(&WHITE);
    let areas = root.split_by_breakpoints([944], [80]);

    let mut scatter_ctx = ChartBuilder::on(&areas[2])
    .x_label_area_size(40)
    .y_label_area_size(40)
    .build_cartesian_2d(-4f64..8f64, -7f64..7f64).unwrap();
    scatter_ctx
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .draw();


    let vertices = vec![
        [0.3, 0.1],
        [1.0, 0.2],
        [0.1, 1.0],
        [0.5, 0.5]
    ];
    let mut delaunay_tree = DelaunayTree::<2,3>::new(vertices.clone());
    assert_eq!(delaunay_tree.max_simplex_id, 3);
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(0.0..1.0);

    // println!("{:?}", delaunay_tree.vertices);
    // scatter_ctx.draw_series(
    //     delaunay_tree.vertices
    //         .iter()
    //         .map(|[x, y]| Circle::new((*x, *y), 2, RED.filled())),
    // );
    // root.present();

    // for i in 0..4 {
    //     println!("{:?}", delaunay_tree.simplices);
    //     let update = TreeUpdate::new(vertices[i], &delaunay_tree);
    //     println!("{:?}", update);
    //     delaunay_tree.insert_point(update);
    //     scatter_ctx.draw_series(
    //         delaunay_tree.vertices
    //             .iter()
    //             .map(|[x, y]| Circle::new((*x, *y), 2, RED.filled())),
    //     );
    //     root.present();
    //     delaunay_tree.check_delaunay();
    // }
    for i in 0..100 {
        let point = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
        let update = TreeUpdate::new(point, &delaunay_tree);
        delaunay_tree.insert_point(update);
        delaunay_tree.check_delaunay();
    }
}
