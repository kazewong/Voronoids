using Revise
using Voronoids
using TimerOutputs
using BenchmarkTools
using ThreadPinning
using Random

pinthreads(:cores)

const tmr = TimerOutput()

Random.seed!(1234)

n = 1000
n_dims = 3

test_points = [rand(n_dims) for i in 1:n]
test_points2 = [rand(n_dims) for i in 1:1e7]
if n_dims == 3
    tree = initialize_tree_3d(test_points)
else
    tree = initialize_tree_2d(test_points)
end
for i in 1:n
    insert_point(tree, test_points[i])
end

test_tree = deepcopy(tree)

new_n = 100

sites, index = identify_nonconflict_points(test_points2[1:new_n], test_tree)
batch_insert_point(test_points2[1:new_n], test_tree)

for i in 1:new_n
    if index[i]
        insert_point(tree, test_points2[i])
    end
end

# check_delaunay(tree)