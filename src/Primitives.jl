using AdaptiveKDTrees.KNN
struct Vertex
    id::Int
    position::Vector{Float64}
end

function initialize_vertex(n::Int; n_dims::Int = 3)::Vector{Vertex}
    return [Vertex(i, rand(n_dims)) for i in 1:n]
end

mutable struct DelaunayTreeNode
    id::Int
    dead::Bool
    vertices::Vector{Int}
    circumcenter::Vector{Float64}
    radius::Float64
end

mutable struct DelaunayTree
    id::Vector{Int}
    simplices::Vector{Vector{Int}}
    dead::Vector{Bool}
    centers::Vector{Vector{Float64}}
    radii::Vector{Float64}
    parent_relation::Vector{Int}
    children_relation::Vector{Vector{Int}}
    step_children_relation::Vector{Dict{Vector{Int},Vector{Int}}}
    neighbors_relation::Vector{Vector{Int}}

    vertices::Vector{Vector{Float64}}

    # simplices::Dict{Int, DelaunayTreeNode}

    # step_children_relation::Dict{Int,Dict{Vector{Int},Vector{Int}}}
end

mutable struct KDDelaunayTree
    vertices::Vector{Vector{Float64}}
    dead::Vector{Bool}
    kdtree::KDTree
    simplices::Vector{Vector{Int}}
    vertices_simplex::Vector{Vector{Int}}
    centers::Vector{Vector{Float64}}
    radii::Vector{Float64}
    neighbors_relation::Vector{Vector{Int}}
end
