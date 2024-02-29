using Revise
using Voronoids
using Plots
using Random

Random.seed!(1234)

n = 10
n_dims = 2

test_points = [rand(n_dims) for i in 1:n]

tree = initialize_tree_2d(test_points)

for i in 1:n
    println(i)
    add_vertex!(tree, test_points[i], n_dims=n_dims)
end

