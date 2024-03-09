using AdaptiveKDTrees.KNN
using BoundingSphere

mutable struct DelaunayTree
    vertices::Vector{Vector{Float64}}
    kdtree::KDTree
    vertices_simplex::Vector{Vector{Int}}

    simplices::Dict{Int, Vector{Int}}
    centers::Dict{Int, Vector{Float64}}
    radii::Dict{Int, Float64}
    neighbors_relation::Dict{Int, Vector{Int}}
    max_simplices_id::Int
end

struct TreeUpdate
    vertices::Vector{Float64} # New vertices
    killed_sites::Vector{Int} # id of the killed simplices
    simplices::Vector{Vector{Int}} # Vertices of the new simplices
    simplices_id::Vector{Int} # Id of the new simplices. Relatively indexed, i.e. max(simplices_id) = length(simplices)
    centers::Vector{Vector{Float64}} # Centers of the circumsphere of the new simplices
    radii::Vector{Float64} # Radii of the circumsphere of the new simplices
    neighbors_id::Vector{Tuple{Int,Int}} # (Neighbor Id, killed site Id)
    new_neighbors_id::Vector{Tuple{Int,Int}} # (New site Id1, New site Id2). Relatively indexed.
end

function insert_point!(tree::DelaunayTree, update::TreeUpdate)
    killed_sites_ids = sort(update.killed_sites)
    push!(tree.vertices, update.vertices)

    # Update neighbor relation with killed sites
    for i in 1:length(update.neighbors_id)
        neighbor_id, killed_id = update.neighbors_id[i]
        tree.neighbors_relation[neighbor_id][tree.neighbors_relation[neighbor_id] .== killed_id] .= tree.max_simplices_id+update.simplices_id[i]
    end

    # Update neighbor relation between new sites
    for i in 1:length(update.simplices_id)
        tree.neighbors_relation[tree.max_simplices_id+update.simplices_id[i]] = [update.neighbors_id[i][1]]
    end
    for i in 1:length(update.new_neighbors_id)
        push!(tree.neighbors_relation[tree.max_simplices_id+update.new_neighbors_id[i][1]], tree.max_simplices_id+update.new_neighbors_id[i][2])
    end

    # Update vertices_simplex
    push!(tree.vertices_simplex, Vector{Int}())
    for i in 1:length(update.simplices)
        for j in update.simplices[i]
            push!(tree.vertices_simplex[j], tree.max_simplices_id+update.simplices_id[i])
        end
    end
    for i in 1:length(killed_sites_ids)
        for j in tree.simplices[killed_sites_ids[i]]
            tree.vertices_simplex[j] = tree.vertices_simplex[j][tree.vertices_simplex[j] .!= killed_sites_ids[i]]
        end
    end

    # Update simplices
    for i in 1:length(update.simplices_id)
        tree.simplices[tree.max_simplices_id+update.simplices_id[i]] = update.simplices[i]
        tree.centers[tree.max_simplices_id+update.simplices_id[i]] = update.centers[i]
        tree.radii[tree.max_simplices_id+update.simplices_id[i]] = update.radii[i]
    end

    # Remove killed sites
    for killed_sites_id in killed_sites_ids
        delete!(tree.simplices, killed_sites_id)
        delete!(tree.centers, killed_sites_id)
        delete!(tree.radii, killed_sites_id)
        delete!(tree.neighbors_relation, killed_sites_id)
    end

    tree.max_simplices_id += maximum(update.simplices_id)
end

function initialize_tree_3d(positions::Vector{Vector{Float64}})::DelaunayTree
    center, radius = boundingsphere(positions)
    radius = radius*5
    first_vertex = center + [0, 0, radius]
    second_vertex = center + [radius * cos(0), radius * sin(0), -radius ]
    third_vertex = center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3), -radius]
    fourth_vertex = center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3), -radius]
    ghost_vertex = [center + [0, 0, radius], center + [0, 0, radius], center + [0, 0, radius], center + [radius  * cos(0), radius * sin(0), -radius]]
    vertices = [first_vertex, second_vertex, third_vertex, fourth_vertex, ghost_vertex...]
    kd_tree = KDTree(reduce(hcat, vertices))
    vertices_simplex = [[1, 2, 3, 4], [1, 2, 4, 5], [1, 2, 3, 5], [1, 3, 4, 5], [2], [3], [4], [5]]

    simplicies = Dict(1 => [1, 2, 3, 4], 2 => [5, 1, 2, 3], 3 => [6, 1, 3, 4], 4 => [7, 1, 4, 2], 5 => [8, 2, 3, 4])
    centers = Dict(1 => center, 2 => [0, 0, 0], 3 => [0, 0, 0], 4 => [0, 0, 0], 5 => [0, 0, 0])
    radii = Dict(1 => radius, 2 => 0, 3 => 0, 4 => 0, 5 => 0)
    neighbors_relation = Dict(1 => [2, 3, 4, 5], 2 => [1], 3 => [1], 4 => [1], 5 => [1])
    

    return DelaunayTree(vertices, kd_tree, vertices_simplex, simplicies, centers, radii, neighbors_relation, 5)
end

function initialize_tree_2d(positions::Vector{Vector{Float64}})::DelaunayTree
    center, radius = boundingsphere(positions)
    radius = radius*3
    first_vertex = center + [radius * cos(0* pi / 3), radius * sin(0 * pi / 3)]
    second_vertex = center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3)]
    third_vertex = center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3)]
    ghost_vertex = [center + [radius * cos(0 * pi / 3), radius * sin(0 * pi / 3)], center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3)], center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3)]]
    vertices = [first_vertex, second_vertex, third_vertex, ghost_vertex...]
    kd_tree = KDTree(reduce(hcat, vertices))
    vertices_simplex = [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2], [3], [4]]

    simplicies = Dict(1 => [1, 2, 3], 2 => [4, 1, 2], 3 => [6, 1, 3], 4 => [5, 2, 3])
    centers = Dict(1 => center, 2 => [0, 0], 3 => [0, 0], 4 => [0, 0])
    radii = Dict(1 => radius, 2 => 0, 3 => 0, 4 => 0)
    neighbors_relation = Dict(1 => [2, 3, 4], 2 => [1], 3 => [1], 4 => [1])

    return DelaunayTree(vertices, kd_tree, vertices_simplex, simplicies, centers, radii, neighbors_relation, 4)
end

export initialize_tree_2d, initialize_tree_3d, DelaunayTree, insert_point!, TreeUpdate