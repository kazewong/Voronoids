using Revise
using Voronoids
using Plots
using Random
using NPZ

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

function make_plot(n_step::Int, tree::DelaunayTree, x_lim::Tuple{Float64, Float64}, y_lim::Tuple{Float64, Float64}; show_circle::Bool=true)
    n = length(tree.vertices[7:7+n_step-1])
    colors = distinguishable_colors(n)
    vertices = tree.vertices[7:7+n_step-1]
    p = scatter(map(x -> x[1], vertices), map(x -> x[2], vertices), c=colors, showaxis=false, legend=false, xlims=x_lim, ylims=y_lim)

    for i in collect(keys(tree.simplices))
        if sum(tree.simplices[i].<6) < 1
            x,y = plot_simplex_2d(i, tree)
            p = plot!(x, y, size=(800, 800), color=:black)
            if show_circle
                x,y = plot_circle(tree.centers[i], tree.radii[i])
                p = plot!(x, y, size=(800, 800), color=:black, alpha=0.3)
            end
        end
    end
    return p
end

function make_frame(n::Int, show_circle::Bool=true; name::String="", plot_center::Vector{Float64}=[0.0, 0.0], plot_radius::Float64=0.0, seed::Int=12345)


    Random.seed!(seed)
    n_dims = 2

    test_points = [rand(n_dims) for i in 1:n]

    tree = initialize_tree_2d(test_points)


    for i in 1:n
        add_vertex!(tree, test_points[i], n_dims=n_dims)
    end

    p = make_plot(n, tree, (plot_center[1]-plot_radius, plot_center[1]+plot_radius), (plot_center[2]-plot_radius, plot_center[2]+plot_radius), show_circle=show_circle)

    if name != ""
        savefig(p, name)
    else
        savefig(p, "frame_$(lpad(n, 4, '0')).png")
    end
end

n_max = 20
seed = 9527

Random.seed!(seed)

n_dims = 2

test_points = [rand(n_dims) for i in 1:n_max]

tree = initialize_tree_2d(test_points)

plot_center = tree.centers[1]
plot_radius = tree.radii[1]*0.5

for i in 1:n_max
    add_vertex!(tree, test_points[i], n_dims=n_dims)
end

npzwrite("tree.npz", Dict("vertices"=>reduce(hcat,tree.vertices), "centers"=>reduce(hcat,collect(values(tree.centers))), "radii"=>collect(values(tree.radii)), "keys"=>collect(keys(tree.simplices)), "values"=>reduce(hcat, collect(values(tree.simplices)))))

# for i in 3:n_max
#     make_frame(i, true, plot_center=plot_center, plot_radius=plot_radius, seed=seed)
# end

# make_frame(n_max, false, name="frame_final_no_circle.png", plot_center=plot_center, plot_radius=plot_radius, seed=seed)