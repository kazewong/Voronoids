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

n_insert = 6553

occupancy = Dict{Int, Vector{Int}}()
lk = ReentrantLock()
channel = Channel{Tuple{Int, Vector{Float64}, Vector{Int}}}(n_parallel)
t1 = Threads.@spawn queue_multiple_points!(channel, test_points2[1:n_insert], occupancy, parallel_tree, lk, batch_size=n_parallel)

t2 = consume_multiple_points!(n_insert, channel, parallel_tree, occupancy, lk, n_dims)
println(length(parallel_tree.vertices))
# t = @async parallel_insert!(test_points2[1:n_parallel], parallel_tree, n_dims=n_dims)

# channel, b, c = fetch(t)

270086
10524
10537