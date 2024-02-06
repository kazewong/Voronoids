struct Vertex
    id::Int
    position::Vector{Float64}
end

function initialize_vertex(n::Int)::Vector{Vertex}
    return [Vertex(i, rand(3)) for i in 1:n]
end

mutable struct DelaunayTreeNode
    id::Int
    dead::Bool
    vertices::Vector{Int}
end

struct DelaunayTree
    vertices::Vector{Vertex}
    simplices::Vector{DelaunayTreeNode}
    children_relation::Dict{Int,Vector{Int}}
    step_children_relation::Dict{Int,Vector{Vector{Int}}}
    neighbors_relation::Dict{Int,Vector{Int}}
end
