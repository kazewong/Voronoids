using Revise
using Voronoids
using TimerOutputs
using BenchmarkTools
using ThreadPinning
using Random

pinthreads(:cores)

const tmr = TimerOutput()

Random.seed!(1234)

n = 3000
n_dims = 3

test_points = [rand(n_dims) for i in 1:n]
test_points2 = [rand(n_dims) for i in 1:1e5]
if n_dims == 3
    tree = initialize_tree_3d(test_points)
else
    tree = initialize_tree_2d(test_points)
end
for i in 1:n
    insert_point(tree, test_points[i])
end
parallel_tree = deepcopy(tree)

function insert_all_points!(tree, points)
    conflicts = Vector{Bool}()
    batch_size = 256
    for i in 1:batch_size:length(test_points2)
        conflicts = vcat(conflicts, batch_insert_point(test_points2[i:min(i+batch_size-1, length(test_points2))], tree))
        println(sum(conflicts))
    end
end

# insert_all_points!(parallel_tree, test_points2)
batch_insert_point(test_points2[1024:2048], parallel_tree)

# insert_point_parallel!(parallel_tree, test_points2, batch_factor=16)
