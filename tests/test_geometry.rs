use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use voronoids::geometry::{bounding_sphere, circumsphere};
#[test]
fn test_circumsphere() {
    let vertices = [
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let (center, radius) = circumsphere(vertices);
    assert_eq!(center, [0.5, 0.5, 0.5]);
    assert_eq!(radius, 0.8660254037844386);
}

#[test]
fn test_boundsphere() {
    const N_TEST: usize = 1000;
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(-1.0..1.0);
    let mut point_test: [[f64; 3]; N_TEST] = [[0.0; 3]; N_TEST];
    for i in 0..N_TEST {
        point_test[i] = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
    }
    let (center, radius) = bounding_sphere(point_test.to_vec());
    println!("Center: {:?}", center);
    println!("Radius: {:?}", radius);
}
