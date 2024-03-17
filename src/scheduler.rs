use rayon::prelude::*;
use std::collections::HashMap;

use crate::delaunay_tree::DelaunayTree;

pub fn make_queue<const N: usize, const M: usize>(
    vertices: Vec<[f64; N]>,
    tree: DelaunayTree<N, M>,
) -> Vec<(usize, [f64; N], Vec<usize>)> {
    vertices.into_par_iter().enumerate().map(|(id, vertex)| {
        let killed_site = tree.locate(vertex);
        let neighbors:Vec<usize> = killed_site.into_iter().map(|site| {
            tree.neighbors.get(&site).unwrap().clone()
        }).flatten().collect();
        (id, vertex, neighbors)
    }).collect()
}
