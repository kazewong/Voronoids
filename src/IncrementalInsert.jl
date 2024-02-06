using BoundingSphere
using LinearAlgebra
include("Primitives.jl")

function initial_bounding_simplex(points::Vector{Vertex})
    positions = map(x -> x.position, points)
    center, radius = boundingsphere(positions)
    first_vertex = Vertex(1, center + [0, 0, radius])
    second_vertex = Vertex(2, center + [radius/sqrt(2) * cos(0), radius/sqrt(2) * sin(0), -radius/sqrt(2)])
    third_vertex = Vertex(3, center + [radius/sqrt(2) * cos(2*pi/3), radius/sqrt(2) * sin(2*pi/3), -radius/sqrt(2)])
    fourth_vertex = Vertex(4, center + [radius/sqrt(2) * cos(4*pi/3), radius/sqrt(2) * sin(4*pi/3), -radius/sqrt(2)])
    return Simplex(1, [first_vertex, second_vertex, third_vertex, fourth_vertex], [], 1, true)
end

function in_sphere(simplex::Simplex, point::Vertex)::Bool
    position = reduce(hcat, map(x -> x.position, simplex.vertices)) .- point.position
    position = vcat(position,mapslices(norm, position, dims=1))
    return det(position) > 0
end

function find_simplex(points::Vertex, simplices::Vector{Simplex})
    
end

