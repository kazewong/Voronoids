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
parallel_tree = deepcopy(tree)

sites = locate(Vector{Int}(), test_points2[1], tree)

for i in 1:n
    add_vertex!(tree, test_points[i], n_dims=n_dims)
end

site_list = identify_conflicts(test_points2[1:1024], tree)