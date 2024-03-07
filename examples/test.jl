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
test_points2 = [rand(n_dims) for i in 1:1e6]
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

n_parallel = 10000

n_insert = 1000000

occupancy = Dict{Int, Vector{Int}}()
queue = make_queue!(test_points2[1:n_insert], occupancy, parallel_tree, batch_size=n_parallel)
placement = find_placement(getindex.(queue, 3), occupancy)
# tconsume_multiple_points!(queue, parallel_tree, occupancy, lk, n_dims)

parallel_insert!(test_points2[1:n_insert], parallel_tree, n_dims=n_dims, batch_size=n_parallel)