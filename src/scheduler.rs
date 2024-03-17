use std::{collections::HashMap, hash::Hash};
use rayon::prelude::*;

use crate::delaunay_tree::DelaunayTree;

pub fn make_queue<const N: usize, const M: usize>(
    vertices: Vec<[f64; N]>,
    tree: DelaunayTree<N, M>,
) -> Vec<(usize, [f64; N], Vec<usize>)> {
    vertices.into_par_iter().enumerate().map(|(id, vertex)| {
        let killed_site = tree.locate(vertex);
        let mut neighbors:Vec<usize> = killed_site.into_iter().map(|site| {
            tree.neighbors.get(&site).unwrap().clone()
        }).flatten().collect();
        neighbors.sort();
        neighbors.dedup();
        (id, vertex, neighbors)
    }).collect()
}

pub fn find_placement(queue: Vec<(usize, [f64; 3], Vec<usize>)>) {//-> Vec<usize>{
    let mut occupancy: HashMap<usize, Vec<usize>> = HashMap::new();
    queue.into_iter().for_each(|(id, vertex, neighbors)| {
        for neighbor in neighbors {
            if occupancy.contains_key(&neighbor) {
                let mut sites = occupancy.get_mut(&neighbor).unwrap();
                sites.push(id);
            } else {
                occupancy.insert(neighbor, vec![id]);
            }
        }
    });
    println!("{:?}", occupancy);
}