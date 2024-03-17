<<<<<<< HEAD
use nalgebra::zero;
=======
>>>>>>> 3a7f414 (Refactor make_queue function in scheduler.rs)
use rayon::prelude::*;
use std::{collections::HashMap, hash::Hash};

use crate::delaunay_tree::DelaunayTree;

pub fn make_queue<const N: usize, const M: usize>(
    vertices: Vec<[f64; N]>,
    tree: &DelaunayTree<N, M>,
) -> Vec<(usize, [f64; N], Vec<usize>)> {
    vertices
        .into_par_iter()
        .enumerate()
        .map(|(id, vertex)| {
            let killed_site = tree.locate(vertex);
            let mut neighbors: Vec<usize> = killed_site
                .into_iter()
                .map(|site| tree.neighbors.get(&site).unwrap().clone())
                .flatten()
                .collect();
<<<<<<< HEAD
            neighbors = neighbors
                .into_iter()
                .map(|site| tree.neighbors.get(&site).unwrap().clone())
                .flatten()
                .collect();
=======
>>>>>>> 3a7f414 (Refactor make_queue function in scheduler.rs)
            neighbors.sort();
            neighbors.dedup();
            (id, vertex, neighbors)
        })
        .collect()
}

<<<<<<< HEAD
pub fn find_placement<const N:usize>(queue: &Vec<(usize, [f64; N], Vec<usize>)>) -> Vec<usize> {
=======
pub fn find_placement(queue: Vec<(usize, [f64; 3], Vec<usize>)>) {
>>>>>>> 3a7f414 (Refactor make_queue function in scheduler.rs)
    //-> Vec<usize>{
    let mut occupancy: HashMap<usize, Vec<usize>> = HashMap::new();
    queue
        .clone()
        .into_iter()
        .for_each(|(id, vertex, neighbors)| {
            for neighbor in neighbors {
                if occupancy.contains_key(&neighbor) {
                    let mut sites = occupancy.get_mut(&neighbor).unwrap();
                    sites.push(id);
                } else {
                    occupancy.insert(neighbor, vec![id]);
                }
            }
        });
    let length = queue.len();
    let mut placement: Vec<usize> = Vec::with_capacity(length);
    for i in 0..length {
        let mut max_distance = 0;
        for neighbor in queue[i].2.clone() {
            let mut current_distance = 0;
            for site in occupancy.get(&neighbor).unwrap() {
                if *site == i {
                    current_distance += 1;
                    break;
                }else{
                    current_distance = std::cmp::max(current_distance, placement[*site]);
                }
            }
            max_distance = std::cmp::max(max_distance, current_distance);
        }
<<<<<<< HEAD
        placement.push(max_distance);
    }
    placement
=======
    });
>>>>>>> 3a7f414 (Refactor make_queue function in scheduler.rs)
}
