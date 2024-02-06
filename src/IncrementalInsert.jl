using Revise
using BoundingSphere
using LinearAlgebra
include("Primitives.jl")

function initialize_tree(points::Vector{Vertex})::DelaunayTree
    positions = map(x -> x.position, points)
    center, radius = boundingsphere(positions)
    first_vertex = Vertex(1, center + [0, 0, radius])
    second_vertex = Vertex(2, center + [radius / sqrt(2) * cos(0), radius / sqrt(2) * sin(0), -radius / sqrt(2)])
    third_vertex = Vertex(3, center + [radius / sqrt(2) * cos(2 * pi / 3), radius / sqrt(2) * sin(2 * pi / 3), -radius / sqrt(2)])
    fourth_vertex = Vertex(4, center + [radius / sqrt(2) * cos(4 * pi / 3), radius / sqrt(2) * sin(4 * pi / 3), -radius / sqrt(2)])
    node = DelaunayTreeNode(1, false, [first_vertex, second_vertex, third_vertex, fourth_vertex])
    return DelaunayTree([node], Dict(1 => []), Dict(1 => []), Dict(1 => []))
end

function in_sphere(node::DelaunayTreeNode, point::Vertex)::Bool
    position = reduce(hcat, map(x -> x.position, node.vertices)) .- point.position
    position = vcat(position, mapslices(norm, position, dims=1))
    return det(position) > 0
end

function locate(visited_ids::Vector{Int}, output::Vector{Int}, vertex::Vertex, current_node::DelaunayTreeNode, tree::DelaunayTree)::Vector{Int}
    if current_node.id ∉ visited_ids && in_sphere(current_node, vertex)
        push!(visited_ids, current_node.id)
        push!(output, current_node.id)
        childrens = tree.nodes[tree.children_relation[current_node.id]]
        check_children= Vector{DelaunayTreeNode}()
        if length(childrens) > 0
            check_children = map(x -> locate(visited_ids, output, vertex, x, tree), childrens)
        end
        step_childrens = tree.step_children_relation[current_node.id]
        check_step_children = []
        if length(step_childrens) > 0
            nodes = tree.nodes[reduce(vcat, step_childrens)]
            check_step_children = map(x -> locate(visited_ids, output, vertex, x, tree), step_childrens)
        end
        return vcat(output, check_children, check_step_children)
    else
        return output
    end
end

# function insert_point(nodes::DelaunayTree, point::Vertex)

# end

test_points = initialize_vertex(100)
root_node = initialize_tree(test_points)
println(locate(Vector{Int}(), Vector{Int}(), test_points[1], root_node.nodes[1], root_node))