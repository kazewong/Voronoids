use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use voronoids::geometry::{bounding_sphere, circumsphere};
#[test]
fn test_circumsphere() {
    let vertices = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let (center, radius) = circumsphere(vertices);
    assert_eq!(center, [0.5, 0.5, 0.5]);
    assert_eq!(radius, 0.8660254037844386);
}

#[test]
fn test_boundsphere() {
    const n_test: usize = 1000;
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(-1.0..1.0);
    let mut point_test: [[f64; 3]; n_test] = [[0.0; 3]; n_test];
    for i in 0..n_test {
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
