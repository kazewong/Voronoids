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
                .map(|site| tree.simplices.get(&site).unwrap().neighbors.clone())
                .flatten()
                .into_iter()
                .map(|site| tree.simplices.get(&site).unwrap().neighbors.clone())
                .flatten()
                .collect();
            neighbors.sort();
            neighbors.dedup();
            (id, vertex, neighbors)
        })
        .collect()
}

pub fn find_placement<const N: usize>(queue: &Vec<(usize, [f64; N], Vec<usize>)>) -> Vec<usize> {
    let mut occupancy: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut placement: Vec<usize> = vec![0; queue.len()];
    queue.into_iter().for_each(|(id, _, neighbors)| {
        for neighbor in neighbors {
            if occupancy.contains_key(&neighbor) {
                let sites = occupancy.get_mut(&neighbor).unwrap();
                sites.push(*id);
            } else {
                occupancy.insert(*neighbor, vec![*id]);
            }
        }
        placement[*id] = neighbors
            .iter()
            .map(|site| {
                if occupancy[site].len() == 1 {
                    1
                } else {
                    placement[occupancy[site][occupancy[site].len() - 2]] + 1
                }
            })
            .max()
            .unwrap()
    });
    placement
}
