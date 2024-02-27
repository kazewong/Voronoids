using AdaptiveKDTrees.KNN
using BoundingSphere

mutable struct DelaunayTree
    vertices::Vector{Vector{Float64}}
    kdtree::KDTree
    vertices_simplex::Vector{Vector{Int}}

    simplices::Vector{Vector{Int}}
    centers::Vector{Vector{Float64}}
    radii::Vector{Float64}
    neighbors_relation::Vector{Vector{Int}}
end

struct TreeUpdate
    vertices::Vector{Float64} # New vertices
    killed_sites::Vector{Int} # id of the killed simplices
    simplices::Vector{Vector{Int}} # Vertices of the new simplices
    centers::Vector{Vector{Float64}} # Centers of the circumsphere of the new simplices
    radii::Vector{Float64} # Radii of the circumsphere of the new simplices
    neighbors_id::Vector{Tuple{Int, Int}} # (Neighbor Id, killed site Id)
    new_neighbors_id::Vector{Tuple{Int, Int}} # (New site Id1, New site Id2)
end

function insert_point!(tree::DelaunayTree, update::TreeUpdate)
    push!(tree.vertices, update.vertices)
    add_point!(tree.kdtree, update.vertices)
    deleteat!(tree.simplices, update.killed_sites)
    deleteat!(tree.centers, update.killed_sites)
    deleteat!(tree.radii, update.killed_sites)
    deleteat!(tree.neighbors_relation, update.killed_sites)
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

    simplicies = [[1, 2, 3, 4], [5, 1, 2, 3], [6, 1, 3, 4], [7, 1, 4, 2], [8, 2, 3, 4]]
    vertices_simplex = [[1,2,3,4], [1, 2, 4, 5], [1, 2, 3, 5], [1, 3, 4, 5], [2], [3], [4], [5]]
    centers = [center, [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]
    radii = [radius, 0, 0, 0, 0]

    neighbors_relation = [[2, 3, 4, 5], [1], [1], [1], [1]]
    kd_tree = KDTree(reduce(hcat, vertices))
    return DelaunayTree(vertices, kd_tree, simplicies, vertices_simplex, centers, radii, neighbors_relation)
end

function initialize_tree_2d(positions::Vector{Vector{Float64}})::DelaunayTree
    center, radius = boundingsphere(positions)
    radius = radius*3
    first_vertex = center + [radius * cos(0* pi / 3), radius * sin(0 * pi / 3)]
    second_vertex = center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3)]
    third_vertex = center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3)]
    ghost_vertex = [center + [radius * cos(0 * pi / 3), radius * sin(0 * pi / 3)], center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3)], center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3)]]
    vertices = [first_vertex, second_vertex, third_vertex, ghost_vertex...]

    simplicies = [[1, 2, 3], [4, 1, 2], [6, 1, 3], [5, 2, 3]]
    vertices_simplex = [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2], [3], [4]]
    centers = [center, [0, 0], [0, 0], [0, 0]]
    radii = [radius, 0, 0, 0]

    neighbors_relation = [[2, 3, 4], [1], [1], [1]]
    kd_tree = KDTree(reduce(hcat, vertices))
    return DelaunayTree(vertices, kd_tree, simplicies, vertices_simplex, centers, radii, neighbors_relation)
end

export initialize_tree_2d, initialize_tree_3d, DelaunayTree