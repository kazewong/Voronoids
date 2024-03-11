use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};

#[test]
fn test_delaunay_tree() {
    let vertices = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let delaunay_tree = DelaunayTree::new(vertices);
    assert_eq!(delaunay_tree.max_simplex_id, 4);
    let new_point = [0.5, 0.5, 0.5];
    let update = TreeUpdate::new(new_point, &delaunay_tree);
    assert!(delaunay_tree.check_delaunay());
}