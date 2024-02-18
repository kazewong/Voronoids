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
for i in 1:100
    insert_point(tree, test_points[i])
end
parallel_tree = deepcopy(tree)

insert_point_parallel!(parallel_tree, test_points2, batch_factor=16)
