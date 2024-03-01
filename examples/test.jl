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
sites = locate(Vector{Int}(), test_points2[1], tree)

for i in 1:n
    add_vertex!(tree, test_points[i], n_dims=n_dims)
end

parallel_tree = deepcopy(tree)

n_parallel = 256


# channel = Channel{Vector{Tuple{Vector{Float64}, Vector{Int}}}}(n_parallel)
# sites, neighbors, occupancy = identify_conflicts(test_points2[1:n_parallel], parallel_tree)
# groups = group_points(site_list, neighbor_list, occupancy)
# t1 = @async queue_multiple_points!(channel, test_points2[1:n_parallel], site_list, groups, parallel_tree, n_dims)
t = @async parallel_insert!(test_points2[1:n_parallel], parallel_tree, n_dims=n_dims)

channel, b, c = fetch(t)

# for i in 1:update_channel.n_avail_items
#     update = take!(update_channel)
#     println("Inserting $(length(update)) vertices")
#     for j in 1:length(update)
#         insert_point!(parallel_tree, update[j])
#     end
#     # insert_point!(tree, update[1])
# end