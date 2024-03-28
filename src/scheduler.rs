use rayon::{prelude::*, vec};
use std::collections::HashMap;

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
                .into_iter()
                .map(|site| tree.neighbors.get(&site).unwrap().clone())
                .flatten()
                .collect();
            neighbors.sort();
            neighbors.dedup();
            (id, vertex, neighbors)
        })
        .collect()
}

pub fn find_placement<const N:usize>(queue: &Vec<(usize, [f64; N], Vec<usize>)>) -> Vec<usize> {
    //-> Vec<usize>{
    let mut occupancy: HashMap<usize, Vec<usize>> = HashMap::new();
    #[cfg(debug_assertions)] println!("Filling occupancy");
    queue.into_iter()
        .for_each(|(id, _, neighbors)| {
            for neighbor in neighbors {
                if occupancy.contains_key(&neighbor) {
                    let sites = occupancy.get_mut(&neighbor).unwrap();
                    sites.push(*id);
                } else {
                    occupancy.insert(*neighbor, vec![*id]);
                }
            }
        });
    let length = queue.len();
    let mut placement: Vec<usize> = vec![0; length];
    #[cfg(debug_assertions)] println!("Finding placement");
    for i in 0..length {
        let mut max_distance = 0;
        queue[i].2.iter().for_each(|neighbor| {
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
        });
        placement[i] = max_distance;
    }
    placement
}
