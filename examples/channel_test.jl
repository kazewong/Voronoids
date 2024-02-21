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
    insert_point(tree, test_points[i])
end

function parallel_locate2(vertices::Vector{Vector{Float64}}, tree::DelaunayTree)::Vector{Vector{Int}}
    output = Vector{Vector{Int}}(undef, length(vertices))
    chunks = collect(Iterators.partition(eachindex(vertices), length(vertices) ÷ Threads.nthreads()))
    Threads.@threads for chunk in chunks
        output[chunk] = batch_locate(vertices[chunk], tree)
    end
    return output
end