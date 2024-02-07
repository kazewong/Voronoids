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
end

struct DelaunayTree
    vertices::Dict{Int,Vertex}
    simplices::Dict{Int, DelaunayTreeNode}
    children_relation::Dict{Int,Vector{Int}}
    step_children_relation::Dict{Int,Dict{Vector{Int},Vector{Int}}}
    neighbors_relation::Dict{Int,Vector{Int}}
end
