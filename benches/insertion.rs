use criterion::{black_box, criterion_group, criterion_main, Criterion};

use rand::{distributions::{Distribution, Uniform}, rngs::StdRng, SeedableRng};
use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};

fn benchmark_locate(c: &mut Criterion){
    const N_POINTS: usize = 10000;
    let mut vertices = vec![];
    let mut rng = StdRng::seed_from_u64(0);
    let dist = Uniform::from(0.0..1.0);
    for _ in 0..N_POINTS {
        let point = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
        vertices.push(point);
    }
    let mut delaunay_tree = DelaunayTree::<3, 4>::new(vertices.clone());
    let new_vertex = [0.5, 0.5, 0.5];
    c.bench_function("locate 10000", |b| b.iter(|| delaunay_tree.locate(new_vertex)));
    c.bench_function("update 10000", |b| b.iter(|| TreeUpdate::new(10001, new_vertex, &delaunay_tree)));
}

criterion_group!(benches, benchmark_locate);
criterion_main!(benches);