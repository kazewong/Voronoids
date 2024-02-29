using Revise
using Voronoids
using Plots
using Random

Random.seed!(1234)

n = 5
n_dims = 2

test_points = [rand(n_dims) for i in 1:n]

tree = initialize_tree_2d(test_points)

for i in 1:n
    add_vertex!(tree, test_points[i], n_dims=n_dims)
end

function plot_simplex_2d(simplex_Id::Int, tree::DelaunayTree)
    x, y = [], []
    vertices = tree.vertices[tree.simplices[simplex_Id]]
    for vertex_id in [1,2,3,1]
        push!(x, vertices[vertex_id][1])
        push!(y, vertices[vertex_id][2])
    end
    return x, y
end

function plot_circle(center::Vector{Float64}, radius::Float64)
    θ = LinRange(0, 2π, 100)
    x = center[1] .+ radius * cos.(θ)
    y = center[2] .+ radius * sin.(θ)
    return x, y
end

scatter(map(x -> x[1], test_points), map(x -> x[2], test_points), c=distinguishable_colors(n)) 

for ids in collect(keys(tree.simplices))
    if all(tree.simplices[ids].>6)
        x,y = plot_simplex_2d(ids, tree)
        p = plot!(x, y, size=(800, 800))
        x,y = plot_circle(tree.centers[ids], tree.radii[ids])
        p = plot!(x, y, size=(800, 800))
    end
end

display(p)
# x,y = plot_simplex_2d(1, tree)
# plot(x, y, label="Points", size=(800, 800))

# p = scatter!(map(x -> x[1], test_points), map(x -> x[2], test_points), label="Points", c=distinguishable_colors(n)) 
# check_delaunay(tree, n_dims=2)