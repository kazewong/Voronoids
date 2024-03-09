use voronoids::geometry::circumsphere;

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