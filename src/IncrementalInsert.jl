using Revise
using BoundingSphere
using LinearAlgebra
using Plots
plotlyjs()
include("Primitives.jl")

function initialize_tree(points::Vector{Vertex})::DelaunayTree
    positions = map(x -> x.position, points)
    center, radius = boundingsphere(positions)
    radius = radius*2
    first_vertex = Vertex(-8, center + [0, 0, radius])
    second_vertex = Vertex(-7, center + [radius * cos(0), radius * sin(0), -radius ])
    third_vertex = Vertex(-6, center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3), -radius])
    fourth_vertex = Vertex(-5, center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3), -radius])
    ghost_vertex = [Vertex(-4, center + [0, 0, radius]), Vertex(-3, center + [0, 0, radius]), Vertex(-2, center + [0, 0, radius]), Vertex(-1, center + [radius  * cos(0), radius * sin(0), -radius])]
    verticies = [first_vertex, second_vertex, third_vertex, fourth_vertex, ghost_vertex...]
    verticies = Dict(map(x -> x.id => x, verticies))
    node = DelaunayTreeNode(1, false, [-8, -7, -6, -5])
    unbounded_node = Vector{DelaunayTreeNode}()
    push!(unbounded_node, DelaunayTreeNode(2, false, [-4, -8, -7, -6]))
    push!(unbounded_node, DelaunayTreeNode(3, false, [-3, -8, -6, -5]))
    push!(unbounded_node, DelaunayTreeNode(4, false, [-2, -8, -5, -7]))
    push!(unbounded_node, DelaunayTreeNode(5, false, [-1, -7, -6, -5]))
    nodes = Dict(1 => node, 2 => unbounded_node[1], 3 => unbounded_node[2], 4 => unbounded_node[3], 5 => unbounded_node[4])
    children_relation = Dict(1 => [], 2 => [], 3 => [], 4 => [], 5 => [])
    step_children_relation = Dict(1 => Dict{Vector{Int},Vector{Int}}(), 2 => Dict{Vector{Int},Vector{Int}}(), 3 => Dict{Vector{Int},Vector{Int}}(), 4 => Dict{Vector{Int},Vector{Int}}(), 5 => Dict{Vector{Int},Vector{Int}}())
    neighbors_relation = Dict(1 => [2, 3, 4, 5], 2 => [1], 3 => [1], 4 => [1], 5 => [1])
    return DelaunayTree(verticies, nodes, children_relation, step_children_relation, neighbors_relation)
end

function in_sphere(node_id::Int, point::Vertex, tree::DelaunayTree)::Bool
    position = reduce(hcat, map(x -> x.position, map(x->tree.vertices[x], tree.simplices[node_id].vertices))) .- point.position
    position = vcat(position, mapslices(norm, position, dims=1))
    return det(position) > 1e-15
end

function locate(visited_ids::Vector{Int}, output::Vector{Int}, vertex::Vertex, current_node_id::Int, tree::DelaunayTree)::Vector{Int}
    if current_node_id ∉ visited_ids && in_sphere(current_node_id, vertex, tree)
        push!(visited_ids, current_node_id)
        # if !tree.simplices[current_node_id].dead
        push!(output, current_node_id)
        # end
        childrens = tree.children_relation[current_node_id]
        for child_id in childrens
            locate(visited_ids, output, vertex, child_id, tree)
        end
        step_childrens = collect(values(tree.step_children_relation[current_node_id]))
        for step_children_id in step_childrens
            locate(visited_ids, output, vertex, step_children_id, tree)
        end
        return output
    else
        return output
    end
end

function common_facet(simplex1::DelaunayTreeNode, simplex2::DelaunayTreeNode)::Vector{Int}
    common = intersect(simplex1.vertices, simplex2.vertices)
    if length(common) == 3
        return common
    else
        return []
    end
end

function insert_point(tree::DelaunayTree, point::Vertex)
    killed_nodes = locate(Vector{Int}(), Vector{Int}(), point, 1, tree)
    println("killed_nodes: ", killed_nodes)
    new_node_id = Vector{Int}()
    for node_id in killed_nodes
        if !tree.simplices[node_id].dead
            tree.simplices[node_id].dead = true
            for neighbor_id in tree.neighbors_relation[node_id]
                if !in_sphere(neighbor_id, point, tree)
                    facet = common_facet(tree.simplices[node_id], tree.simplices[neighbor_id])
                    if length(facet) == 3
                        # Creating new node
                        new_id = length(tree.simplices) + 1
                        println("new_id: ", new_id)
                        println("facet: ", facet)
                        new_node = DelaunayTreeNode(new_id, false, [point.id, facet...])
                        tree.simplices[new_id] = new_node
                        push!(new_node_id, new_node.id)
                        tree.children_relation[new_node.id] = Vector{Int}()
                        tree.step_children_relation[new_node.id] = Dict{Vector{Int}, Vector{Int}}()
                        tree.neighbors_relation[new_node.id] = [neighbor_id]

                        # Updating parent relationship
                        push!(tree.children_relation[node_id],new_node.id)
                        if haskey(tree.step_children_relation[neighbor_id], facet)
                            push!(tree.step_children_relation[neighbor_id][facet], new_node.id)
                        else
                            tree.step_children_relation[neighbor_id][facet] = [new_node.id]
                        end

                        # Updating neighbor relationship for the neighbor of the killed node
                        killed_node_id = findfirst(x->x==node_id, tree.neighbors_relation[neighbor_id])
                        if killed_node_id !== nothing
                            tree.neighbors_relation[neighbor_id][findfirst(x->x==node_id, tree.neighbors_relation[neighbor_id])] = new_node.id
                        end
                    end
                end
            end
        end
    end
    for i in 1:length(new_node_id)
        for j in i+1:length(new_node_id)
            new_id1 = new_node_id[i]
            new_id2 = new_node_id[j]
            facet = common_facet(tree.simplices[new_id1], tree.simplices[new_id2])
            if length(facet) == 3
                push!(tree.neighbors_relation[new_id1], new_id2)
                push!(tree.neighbors_relation[new_id2], new_id1)
            end
        end
    end
    tree.vertices[point.id] = point
end

test_points = initialize_vertex(100)
tree = initialize_tree(test_points)
# println(locate(Vector{Int}(), Vector{Int}(), test_points[1], 1, delaunay_tree))
# insert_point(tree, test_points[1])
# insert_point(tree, test_points[2])

for i in 1:10
    insert_point(tree, test_points[i])
end
# new_tree = deepcopy(tree)


x,y,z = [], [], []
for simplex in values(tree.simplices)
    for vertex_id in simplex.vertices
        push!(x, tree.vertices[vertex_id].position[1])
        push!(y, tree.vertices[vertex_id].position[2])
        push!(z, tree.vertices[vertex_id].position[3])
    end
end

path3d(x, y, z, label="Points", size=(800, 800))
scatter!([test_points[2].position[1]], [test_points[2].position[2]], [test_points[2].position[3]], label="New Point", color="red")