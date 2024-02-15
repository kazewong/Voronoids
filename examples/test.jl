using Voronoids
using TimerOutputs
using BenchmarkTools
using ThreadPinning

pinthreads(:cores)

const tmr = TimerOutput()

n = 10000

test_points = [rand(3) for i in 1:n]
test_points2 = [rand(3) for i in 1:n]
tree = initialize_tree_3d(test_points)
for i in 1:n
    insert_point(tree, test_points[i])
end

# check_delaunay(tree)