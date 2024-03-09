use voronoids::delaunay_tree::DelaunayTree;



#[test]
fn test_delaunay_tree() {
    let vertices = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let delaunay_tree = DelaunayTree::new(vertices);
    assert_eq!(delaunay_tree.max_simplex_id, 5);
    let a = delaunay_tree.locate([0.5, 0.5, 0.5]);
    assert!(delaunay_tree.check_delaunay());
}