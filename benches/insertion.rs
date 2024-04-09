use criterion::{criterion_group, criterion_main, Criterion};

use rand::{
    distributions::{Distribution, Uniform},
    rngs::StdRng,
    SeedableRng,
};
use voronoids::delaunay_tree::{DelaunayTree, TreeUpdate};
use voronoids::scheduler::make_queue;

fn benchmark_locate(c: &mut Criterion) {
    const N_POINTS: usize = 10;
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
    let n_points = delaunay_tree.vertices.len();
    for i in 0..N_POINTS {
        let update = TreeUpdate::new(n_points + i, vertices[i], &delaunay_tree);
        delaunay_tree.insert_point(&update);
    }
    let new_vertex = [0.5, 0.5, 0.5];
    const N_TEST_POINTS: usize = 10000;
    let mut vertices2: Vec<[f64; 3]> = vec![];
    for _ in 0..N_TEST_POINTS {
        let point = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
        vertices2.push(point);
    }
    let queue = make_queue(vertices2.clone(), &delaunay_tree);
    let mut group = c.benchmark_group("insertion_group");
    group.significance_level(0.1).sample_size(10);

    // c.bench_function("locate 10000", |b| b.iter(|| delaunay_tree.locate(new_vertex)));
    // c.bench_function("update 10000", |b| b.iter(|| TreeUpdate::new(10001, new_vertex, &delaunay_tree)));
    // c.bench_function("make_queue 10000", |b| b.iter(|| make_queue(vertices2.clone(), &delaunay_tree)));
    group.bench_function("find_placement 10000", |b| {
        b.iter(|| voronoids::scheduler::find_placement(&queue))
    });
}

fn benchmark_geometry(c: &mut Criterion) {
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
    let n_points = delaunay_tree.vertices.len();
    for i in 0..N_POINTS {
        let update = TreeUpdate::new(n_points + i, vertices[i], &delaunay_tree);
        delaunay_tree.insert_point(&update);
    }
    let new_vertex = [0.5, 0.5, 0.5];
    const N_TEST_POINTS: usize = 10000;
    let mut vertices2: Vec<[f64; 3]> = vec![];
    for _ in 0..N_TEST_POINTS {
        let point = [
            dist.sample(&mut rng),
            dist.sample(&mut rng),
            dist.sample(&mut rng),
        ];
        vertices2.push(point);
    }
    let queue = make_queue(vertices2.clone(), &delaunay_tree);
    c.bench_function("locate 10000", |b| {
        b.iter(|| delaunay_tree.locate(new_vertex))
    });
    c.bench_function("update 10000", |b| {
        b.iter(|| TreeUpdate::new(10001, new_vertex, &delaunay_tree))
    });
    c.bench_function("make_queue 10000", |b| {
        b.iter(|| make_queue(vertices2.clone(), &delaunay_tree))
    });
    c.bench_function("find_placement 10000", |b| {
        b.iter(|| voronoids::scheduler::find_placement(&queue))
    });
}

criterion_group!(benches, benchmark_locate);
criterion_main!(benches);
