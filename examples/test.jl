using Revise
using Voronoids
using TimerOutputs
using BenchmarkTools
using ThreadPinning
using Random

pinthreads(:cores)

const tmr = TimerOutput()

Random.seed!(1234)

n = 10000
n_dims = 3

test_points = [rand(n_dims) for i in 1:n]
test_points2 = [rand(n_dims) for i in 1:1e5]
if n_dims == 3
    tree = initialize_tree_3d(test_points)
else
    tree = initialize_tree_2d(test_points)
end

for i in 1:n
    add_vertex!(tree, test_points[i], n_dims=n_dims)
end

parallel_insert!(test_points2, tree, n_dims=n_dims, batch_size=10000)

parallel_tree = deepcopy(tree)

test_points3 = [rand(n_dims) for i in 1:1e5]
occupancy = Dict{Int, Vector{Int}}()
queue = make_queue!(test_points3, occupancy, parallel_tree)
placement = find_placement(getindex.(queue, 3), occupancy)
all_vertices = getindex.(queue, 2)
test_vertices = all_vertices[findall(x->x==1, placement)]
updates = Vector{TreeUpdate}(undef, length(test_vertices))
n_vertex = length(parallel_tree.vertices)
max_time = 0.0
min_time = Inf
for i in 1:length(test_vertices)
    timer = time()
    updates[i] = make_update(i+n_vertex, test_vertices[i], parallel_tree, n_dims=n_dims)
    max_time = max(max_time, time()-timer)
    min_time = min(min_time, time()-timer)
end

partition = collect(Iterators.partition(1:length(test_vertices), max(length(test_vertices) ÷ Threads.nthreads(),1)))
Threads.@threads for range in partition
    for i in range
        updates[i] = make_update(i+n_vertex, test_vertices[i], parallel_tree, n_dims=n_dims)
    end
end