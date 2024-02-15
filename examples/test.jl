using Voronoids
using TimerOutputs
using BenchmarkTools
using ThreadPinning

pinthreads(:cores)

const tmr = TimerOutput()

n = 10000
n_dims = 3

test_points = [rand(n_dims) for i in 1:n]
test_points2 = [rand(n_dims) for i in 1:n]
if n_dims == 3
    tree = initialize_tree_3d(test_points)
else
    tree = initialize_tree_2d(test_points)
end
for i in 1:n
    insert_point(tree, test_points[i])
end

# check_delaunay(tree)