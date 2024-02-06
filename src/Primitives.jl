struct Vertex
    id::Int
    position::Vector{Float64}
end

function initialize_vertex(n::Int)::Vector{Vertex}
    return [Vertex(i, rand(3)) for i in 1:n]
end

struct Simplex
    vertices::Vector{Vertex}
end

function across(simplex::Simplex)::Simplex
end

function rotate_clockwise(simplex::Simplex)::Simplex
end

struct DelaunayTreeNode
    id::Int
    dead::Bool
    vertices::Vector{Vertex}
end

struct DelaunayTree
    nodes::Vector{DelaunayTreeNode}
    children_relation::Dict{Int, Vector{Int}}
    step_children_relation::Dict{Int, Vector{Vector{Int}}}
    neighbors_relation::Dict{Int, Vector{Int}}
end
