struct Vertex
    id::Int
    position::Vector{Float64}
end

function initialize_vertex(n::Int)::Vector{Vertex}
    return [Vertex(i, rand(3)) for i in 1:n]
end

struct Simplex
    id::Int
    vertices::Vector{Vertex}
    neighbors::Vector{Int}
    orientation::Int
end

function across(simplex::Simplex)::Simplex
end

function rotate_clockwise(simplex::Simplex)::Simplex
end


struct DelaunayTree
    children::Vector{DelaunayTree}
    step_children::Vector{Vector{DelaunayTree}}
    dead::Bool
    simplex::Simplex
end
