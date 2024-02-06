using BoundingSphere
using LinearAlgebra
include("Primitives.jl")

function initial_bounding_simplex(points::Vector{Vertex})::DelaunayTree
    positions = map(x -> x.position, points)
    center, radius = boundingsphere(positions)
    first_vertex = Vertex(1, center + [0, 0, radius])
    second_vertex = Vertex(2, center + [radius/sqrt(2) * cos(0), radius/sqrt(2) * sin(0), -radius/sqrt(2)])
    third_vertex = Vertex(3, center + [radius/sqrt(2) * cos(2*pi/3), radius/sqrt(2) * sin(2*pi/3), -radius/sqrt(2)])
    fourth_vertex = Vertex(4, center + [radius/sqrt(2) * cos(4*pi/3), radius/sqrt(2) * sin(4*pi/3), -radius/sqrt(2)])
    simplex = Simplex(1, [first_vertex, second_vertex, third_vertex, fourth_vertex], [], 1)
    return DelaunayTree([], [], false, simplex)
end

function in_sphere(simplex::Simplex, point::Vertex)::Bool
    position = reduce(hcat, map(x -> x.position, simplex.vertices)) .- point.position
    position = vcat(position,mapslices(norm, position, dims=1))
    return det(position) > 0
end

function locate!(output::Vector{DelaunayTree}, vertex::Vertex, current_node::DelaunayTree)::Vector{DelaunayTree}
    if in_sphere(current_node.simplex, vertex)
        if !current_node.dead        
            current_node.dead = true
            push!(output, current_node)
        end
        list_of_children = map(x -> locate!(output, vertex, x), current_node.children)
        list_of_stepchildren = map(x -> locate!(output, vertex, x), current_node.step_children)
        return vcat(list_of_children, list_of_stepchildren)
    else
        return output
    end
end

function replace_boundary(simplex1::Simplex, simplex2::Simplex, vertex::Vertex)

end

function insert_point(simplices::Vector{Simplex}, point::Vertex)

end

test_points = initialize_vertex(100)
bounding_simplex = initial_bounding_simplex(test_points)