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
    ghost_vertex = [Vertex(5, center + [0, 0, radius]), Vertex(6, center + [0, 0, radius]), Vertex(7, center + [0, 0, radius]), Vertex(8, center + [radius / sqrt(2) * cos(0), radius / sqrt(2) * sin(0), -radius / sqrt(2)])]
    verticies = [first_vertex, second_vertex, third_vertex, fourth_vertex, ghost_vertex...]
    node = DelaunayTreeNode(1, false, [1, 2, 3, 4])
    unbounded_node = Vector{DelaunayTreeNode}()
    push!(unbounded_node, DelaunayTreeNode(2, false, [5, 1, 2, 3]))
    push!(unbounded_node, DelaunayTreeNode(3, false, [6, 1, 2, 4]))
    push!(unbounded_node, DelaunayTreeNode(4, false, [7, 1, 3, 4]))
    push!(unbounded_node, DelaunayTreeNode(5, false, [8, 2, 3, 4]))
    children_relation = Dict(1 => [], 2 => [], 3 => [], 4 => [], 5 => [])
    step_children_relation = Dict(1 => [], 2 => [], 3 => [], 4 => [], 5 => [])
    neighbors_relation = Dict(1 => [2, 3, 4, 5], 2 => [1], 3 => [1], 4 => [1], 5 => [1])
    return DelaunayTree(verticies, [node, unbounded_node...], children_relation, step_children_relation, neighbors_relation)
end

function in_sphere(node_id::Int, point::Vertex, tree::DelaunayTree)::Bool
    position = reduce(hcat, map(x -> x.position, tree.vertices[tree.simplices[node_id].vertices])) .- point.position
    position = vcat(position, mapslices(norm, position, dims=1))
    return abs(det(position)) > 1e-15
end

function locate(visited_ids::Vector{Int}, output::Vector{Int}, vertex::Vertex, current_node_id::Int, tree::DelaunayTree)::Vector{Int}
    if current_node_id ∉ visited_ids && in_sphere(current_node_id, vertex, tree)
        push!(visited_ids, current_node_id)
        if !tree.simplices[current_node_id].dead
            push!(output, current_node_id)
        end
        childrens = tree.simplices[tree.children_relation[current_node_id]]
        check_children = Vector{DelaunayTreeNode}()
        if length(childrens) > 0
            check_children = map(x -> locate(visited_ids, output, vertex, x, tree), childrens)
        end
        step_childrens = tree.step_children_relation[current_node_id]
        check_step_children = []
        if length(step_childrens) > 0
            nodes = tree.simplices[reduce(vcat, step_childrens)]
            check_step_children = map(x -> locate(visited_ids, output, vertex, x, tree), step_childrens)
        end
        return vcat(output, check_children, check_step_children)
    else
        return output
    end
end

function insert_point(tree::DelaunayTree, point::Vertex)
    killed_nodes = locate(Vector{Int}(), Vector{Int}(), point, tree.simplices[1], tree)
    for node_id in killed_nodes
        tree.simplices[node_id].dead = true
        for neighbor_id in tree.neighbors_relation[node_id]
            if !in_sphere(tree.simplices[neighbor_id], point)
                continue
            end
        end
    end
end

test_points = initialize_vertex(100)
delaunay_tree = initialize_tree(test_points)
println(locate(Vector{Int}(), Vector{Int}(), test_points[1], 1, delaunay_tree))