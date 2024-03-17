use plotters::prelude::*;
use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};

const OUT_FILE_NAME: &str = "./test.png";

#[test]
fn test_delaunay_tree3D() {
    let mut vertices = vec![];
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(0.0..1.0);
    for _ in 0..10 {
        let point = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
        vertices.push(point);
    }
    let mut delaunay_tree = DelaunayTree::<3, 4>::new(vertices.clone());
    assert_eq!(delaunay_tree.max_simplex_id, 4);


    for i in 0..5 {
        let update = TreeUpdate::new(vertices[i], &delaunay_tree);
        delaunay_tree.insert_point(update);
        delaunay_tree.check_delaunay();
        println!("{:?}", delaunay_tree.vertices);
        let mut keys = delaunay_tree.simplices.keys().collect::<Vec<_>>();
        keys.sort();
        for key in keys {
            println!(
                "{:?}, {:?}, {:?}",
                &key, delaunay_tree.simplices[&key], delaunay_tree.neighbors[&key]
            );
        }

    }
    let root: DrawingArea<BitMapBackend<'_>, plotters::coord::Shift> = BitMapBackend::gif("test.gif", (1024, 768),500).unwrap().into_drawing_area();
    for i in 1..10 {
        root.fill(&WHITE);
        let mut chart = ChartBuilder::on(&root)
            .caption("3D Plot Test".to_string(), ("sans", 20))
            .build_cartesian_3d(-1.0..1.0,-1.0..1.0, -1.0..1.0).unwrap();
        chart
            .configure_axes()
            .light_grid_style(BLACK.mix(0.15))
            .max_light_lines(3)
            .draw();
        chart.with_projection(|mut p| {
            p.pitch =   0.0* (2.0*std::f64::consts::PI); // 90 degreen pitch, thus we are looking the plot from top
            p.yaw = (i as f64/40.0)*(2.0*std::f64::consts::PI);      // Make plot's X axis parallel to screen's X axis    
            p.into_matrix() // build the projection matrix
        });
        for (key, value) in delaunay_tree.simplices.iter() {
            if value.iter().all(|x| *x > 7) {
                let _ = chart.draw_series(std::iter::once(PathElement::new(
                    vec![
                        (
                            delaunay_tree.vertices[value[0]][0],
                            delaunay_tree.vertices[value[0]][1],
                            delaunay_tree.vertices[value[0]][2],
                        ),
                        (
                            delaunay_tree.vertices[value[1]][0],
                            delaunay_tree.vertices[value[1]][1],
                            delaunay_tree.vertices[value[1]][2],
                        ),
                        (
                            delaunay_tree.vertices[value[2]][0],
                            delaunay_tree.vertices[value[2]][1],
                            delaunay_tree.vertices[value[2]][2],
                        ),
                        (
                            delaunay_tree.vertices[value[3]][0],
                            delaunay_tree.vertices[value[3]][1],
                            delaunay_tree.vertices[value[3]][2],
                        ),
                        (
                            delaunay_tree.vertices[value[0]][0],
                            delaunay_tree.vertices[value[0]][1],
                            delaunay_tree.vertices[value[0]][2],
                        ),
                        (
                            delaunay_tree.vertices[value[2]][0],
                            delaunay_tree.vertices[value[2]][1],
                            delaunay_tree.vertices[value[2]][2],
                        ),
                        (
                            delaunay_tree.vertices[value[1]][0],
                            delaunay_tree.vertices[value[1]][1],
                            delaunay_tree.vertices[value[1]][2],
                        ),
                        (
                            delaunay_tree.vertices[value[3]][0],
                            delaunay_tree.vertices[value[3]][1],
                            delaunay_tree.vertices[value[3]][2],
                        ),
                    ],
                    ShapeStyle {
                        color: BLUE.mix(0.5),
                        filled: true,
                        stroke_width: 1,
                    },
                )));
            }
        }
        root.present();
    }

}


#[test]
fn test_delaunay_tree2D() {
    let root = BitMapBackend::new(OUT_FILE_NAME, (1024, 768)).into_drawing_area();

    root.fill(&WHITE);
    let areas = root.split_by_breakpoints([944], [80]);

    let mut scatter_ctx = ChartBuilder::on(&areas[2])
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-0.2f64..1.2f64, -0.2f64..1.2f64)
        .unwrap();
    let _ = scatter_ctx
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .draw();

    let vertices = vec![[0.3, 0.1], [1.0, 0.2], [0.1, 1.0], [0.5, 0.5]];
    let mut delaunay_tree = DelaunayTree::<2, 3>::new(vertices.clone());
    assert_eq!(delaunay_tree.max_simplex_id, 3);
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(0.0..1.0);

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
        let point = [dist.sample(&mut rng), dist.sample(&mut rng)];
        let update = TreeUpdate::new(point, &delaunay_tree);
        delaunay_tree.insert_point(update);
    }

    delaunay_tree.check_delaunay();

    let _ = scatter_ctx.draw_series(
        delaunay_tree
            .vertices
            .iter()
            .map(|[x, y]| Circle::new((*x, *y), 2, RED.filled())),
    );
    for (key, value) in delaunay_tree.simplices.iter() {
        if value.iter().all(|x| *x > 5) {
            let _ = scatter_ctx.draw_series(std::iter::once(PathElement::new(
                vec![
                    (
                        delaunay_tree.vertices[value[0]][0],
                        delaunay_tree.vertices[value[0]][1],
                    ),
                    (
                        delaunay_tree.vertices[value[1]][0],
                        delaunay_tree.vertices[value[1]][1],
                    ),
                    (
                        delaunay_tree.vertices[value[2]][0],
                        delaunay_tree.vertices[value[2]][1],
                    ),
                    (
                        delaunay_tree.vertices[value[0]][0],
                        delaunay_tree.vertices[value[0]][1],
                    ),
                ],
                ShapeStyle {
                    color: BLUE.mix(0.5),
                    filled: true,
                    stroke_width: 1,
                },
            )));
            // scatter_ctx.draw_series(std::iter::once(Circle::new(
            //     (delaunay_tree.centers[key][0], delaunay_tree.centers[key][1]),
            //     delaunay_tree.radii[key],
            //     ShapeStyle {
            //         color: RED.mix(0.5),
            //         filled: true,
            //         stroke_width: 1,
            //     },
            // )));
        }
    }

}
