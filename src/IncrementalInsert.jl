using Revise
using BoundingSphere
using LinearAlgebra
using Plots
plotlyjs()
using Random
include("Primitives.jl")

function plot_simplex_2d(simplex::DelaunayTreeNode, vertices::Dict{Int, Vertex})
    x, y = [], []
    for vertex_id in [1,2,3,1]
        push!(x, vertices[simplex.vertices[vertex_id]].position[1])
        push!(y, vertices[simplex.vertices[vertex_id]].position[2])
    end
    return x, y
end

function plot_simplex_3d(simplex::DelaunayTreeNode, vertices::Dict{Int, Vertex})
    x, y, z = [], [], []
    for vertex_id in [1,2,3,1,4,3,2,4]
        push!(x, vertices[simplex.vertices[vertex_id]].position[1])
        push!(y, vertices[simplex.vertices[vertex_id]].position[2])
        push!(z, vertices[simplex.vertices[vertex_id]].position[3])
    end
    return x, y, z
end

function sphere(C, r)   # r: radius; C: center [cx,cy,cz]
    n = 100
    u = range(-π, π; length = n)
    v = range(0, π; length = n)
    x = C[1] .+ r*cos.(u) * sin.(v)'
    y = C[2] .+ r*sin.(u) * sin.(v)'
    z = C[3] .+ r*ones(n) * cos.(v)'
    return x, y, z
end

function circumcircle(node_id::Int, tree::DelaunayTree)
    vertices = map(x->tree.vertices[x].position, tree.simplices[node_id].vertices)
    x1, y1 = vertices[1]
    x2, y2 = vertices[2]
    x3, y3 = vertices[3]

    # Midpoints of AB and BC
    D = ((x1 + x2) / 2, (y1 + y2) / 2)
    E = ((x2 + x3) / 2, (y2 + y3) / 2)

    # Slopes of AB and BC
    mAB = (y2 - y1) / (x2 - x1)
    mBC = (y3 - y2) / (x3 - x2)

    # Slopes of perpendicular bisectors
    mD = -1 / mAB
    mE = -1 / mBC

    # Calculating the circumcenter (X, Y)
    X = (mD * D[1] - mE * E[1] + E[2] - D[2]) / (mD - mE)
    Y = mD * (X - D[1]) + D[2]

    # Radius of the circumcircle
    R = sqrt((X - x1)^2 + (Y - y1)^2)

    return ((X, Y), R)
end

function circumsphere(node_id::Int, tree::DelaunayTree)
    vertices = map(x->tree.vertices[x].position, tree.simplices[node_id].vertices)
    # Convert vertices into Julia tuples if they're not already
    v1 = vertices[1]
    v2 = vertices[2]
    v3 = vertices[3]
    v4 = vertices[4]

    if (v1==v2) || (v1==v3) || (v1==v4) || (v2==v3) || (v2==v4) || (v3==v4)
        return ((0, 0, 0), 0)
    end

    a = det([v1[1] v1[2] v1[3] 1;
         v2[1] v2[2] v2[3] 1;
         v3[1] v3[2] v3[3] 1;
         v4[1] v4[2] v4[3] 1])

    Dx = det([v1[1]^2 + v1[2]^2 + v1[3]^2 v1[2] v1[3] 1;
        v2[1]^2 + v2[2]^2 + v2[3]^2 v2[2] v2[3] 1;
        v3[1]^2 + v3[2]^2 + v3[3]^2 v3[2] v3[3] 1;
        v4[1]^2 + v4[2]^2 + v4[3]^2 v4[2] v4[3] 1])

    Dy = - det([v1[1]^2 + v1[2]^2 + v1[3]^2 v1[1] v1[3] 1;
        v2[1]^2 + v2[2]^2 + v2[3]^2 v2[1] v2[3] 1;
        v3[1]^2 + v3[2]^2 + v3[3]^2 v3[1] v3[3] 1;
        v4[1]^2 + v4[2]^2 + v4[3]^2 v4[1] v4[3] 1])

    Dz = det([v1[1]^2 + v1[2]^2 + v1[3]^2 v1[1] v1[2] 1;
        v2[1]^2 + v2[2]^2 + v2[3]^2 v2[1] v2[2] 1;
        v3[1]^2 + v3[2]^2 + v3[3]^2 v3[1] v3[2] 1;
        v4[1]^2 + v4[2]^2 + v4[3]^2 v4[1] v4[2] 1])

    c = det([v1[1]^2 + v1[2]^2 + v1[3]^2 v1[1] v1[2] v1[3];
        v2[1]^2 + v2[2]^2 + v2[3]^2 v2[1] v2[2] v2[3];
        v3[1]^2 + v3[2]^2 + v3[3]^2 v3[1] v3[2] v3[3];
        v4[1]^2 + v4[2]^2 + v4[3]^2 v4[1] v4[2] v4[3]])

    radius = sqrt(Dx^2 + Dy^2 + Dz^2 - 4*a*c) / (2*abs(a))

    return ([Dx/2/a,Dy/2/a,Dz/2/a], radius) # Return the center coordinates and the radius
end


function initialize_tree_3d(points::Vector{Vertex})::DelaunayTree
    positions = map(x -> x.position, points)
    center, radius = boundingsphere(positions)
    radius = radius*3
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

function initialize_tree_2d(points::Vector{Vertex})::DelaunayTree
    positions = map(x -> x.position, points)
    center, radius = boundingsphere(positions)
    radius = radius*3
    first_vertex = Vertex(-6, center + [radius * cos(0* pi / 3), radius * sin(0 * pi / 3)])
    second_vertex = Vertex(-5, center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3)])
    third_vertex = Vertex(-4, center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3)])
    ghost_vertex = [Vertex(-3, center + [radius * cos(0 * pi / 3), radius * sin(0 * pi / 3)]), Vertex(-2, center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3)]), Vertex(-1, center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3)])]
    verticies = [first_vertex, second_vertex, third_vertex, ghost_vertex...]
    verticies = Dict(map(x -> x.id => x, verticies))
    node = DelaunayTreeNode(1, false, [-6, -5, -4])
    unbounded_node = Vector{DelaunayTreeNode}()
    push!(unbounded_node, DelaunayTreeNode(2, false, [-3, -6, -5]))
    push!(unbounded_node, DelaunayTreeNode(3, false, [-1, -6, -4]))
    push!(unbounded_node, DelaunayTreeNode(4, false, [-2, -5, -4]))
    nodes = Dict(1 => node, 2 => unbounded_node[1], 3 => unbounded_node[2], 4 => unbounded_node[3])
    children_relation = Dict(1 => [], 2 => [], 3 => [], 4 => [])
    step_children_relation = Dict(1 => Dict{Vector{Int},Vector{Int}}(), 2 => Dict{Vector{Int},Vector{Int}}(), 3 => Dict{Vector{Int},Vector{Int}}(), 4 => Dict{Vector{Int},Vector{Int}}())
    neighbors_relation = Dict(1 => [2, 3, 4], 2 => [1], 3 => [1], 4 => [1])
    return DelaunayTree(verticies, nodes, children_relation, step_children_relation, neighbors_relation)
end

function in_sphere(node_id::Int, point::Vertex, tree::DelaunayTree; n_dims::Int=3)::Bool
    if n_dims==3
        center, radius = circumsphere(node_id, tree)
        return norm(point.position .- center) < radius
    elseif n_dims==2
        center, radius = circumcircle(node_id, tree)
        return norm(point.position .- center) < radius
    end
end

function locate(visited_ids::Vector{Int}, output::Vector{Int}, vertex::Vertex, current_node_id::Int, tree::DelaunayTree; n_dims::Int = 3)::Vector{Int}
    if current_node_id ∉ visited_ids && in_sphere(current_node_id, vertex, tree, n_dims=n_dims)
        push!(visited_ids, current_node_id)
        push!(output, current_node_id)
        childrens = tree.children_relation[current_node_id]
        for child_id in childrens
            locate(visited_ids, output, vertex, child_id, tree, n_dims=n_dims)
        end
        step_childrens = collect(values(tree.step_children_relation[current_node_id]))
        if length(step_childrens) > 0
            step_childrens = vcat(step_childrens...)
        end
        for step_children_id in step_childrens
            locate(visited_ids, output, vertex, step_children_id, tree, n_dims=n_dims)
        end
        return output
    else
        return output
    end
end

function common_facet(simplex1::DelaunayTreeNode, simplex2::DelaunayTreeNode; n_dims::Int = 3)::Vector{Int}
    common = intersect(simplex1.vertices, simplex2.vertices)
    if length(common) == n_dims
        return common
    else
        return []
    end
end

function insert_point(tree::DelaunayTree, point::Vertex; n_dims::Int=3)
    killed_nodes = locate(Vector{Int}(), Vector{Int}(), point, 1, tree, n_dims=n_dims)
    # println("killed_nodes: ", killed_nodes)
    new_node_id = Vector{Int}()
    for node_id in killed_nodes
        if !tree.simplices[node_id].dead
            tree.simplices[node_id].dead = true
            for neighbor_id in tree.neighbors_relation[node_id]
                # println(in_sphere(neighbor_id, point, tree, n_dims=n_dims))
                if !in_sphere(neighbor_id, point, tree, n_dims=n_dims)
                    # println("neighbor_id: ", neighbor_id)
                    facet = common_facet(tree.simplices[node_id], tree.simplices[neighbor_id], n_dims=n_dims)
                    # println("facet: ", facet)
                    if length(facet) == n_dims
                        # Creating new node
                        new_id = length(tree.simplices) + 1
                        # println("new_id: ", new_id)
                        # println("facet: ", facet)
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
            facet = common_facet(tree.simplices[new_id1], tree.simplices[new_id2], n_dims=n_dims)
            if length(facet) == n_dims
                push!(tree.neighbors_relation[new_id1], new_id2)
                push!(tree.neighbors_relation[new_id2], new_id1)
            end
        end
    end
    tree.vertices[point.id] = point
end

function check_delaunay(tree::DelaunayTree; n_dims::Int=3)
    for i in keys(tree.simplices)
        if !tree.simplices[i].dead
            for j in keys(tree.vertices)
                if j > 0
                    if in_sphere(i, tree.vertices[j], tree, n_dims=n_dims) && j ∉ tree.simplices[i].vertices && all(tree.simplices[i].vertices.>0)
                        println("Error, point ", j, " is inside the circumcircle of simplex ", i)
                    end
                end
            end
        end
    end
end

function test_2d(n::Int; seed::Int)
    Random.seed!(seed)
    n_dims = 2
    test_points = initialize_vertex(n, n_dims=n_dims)
    tree = initialize_tree_2d(test_points)
    
    for point in test_points
        insert_point(tree, point, n_dims=n_dims)
    end
    
    x,y = plot_simplex_2d(tree.simplices[1], tree.vertices)
    p = plot(x, y, label="Points", size=(800, 800))
    for i in 2:length(tree.simplices)
        if !tree.simplices[i].dead && all(tree.simplices[i].vertices.>0)
            x,y = plot_simplex_2d(tree.simplices[i], tree.vertices)
            plot!(x, y, label="Points", size=(800, 800))
        end
    end
    
    scatter!([x for x in map(x -> x.position[1], test_points)], [y for y in map(x -> x.position[2], test_points)], label="Points", color=["red","blue","green", "yellow", "black"])
    
    check_delaunay(tree, n_dims=2)

    return p
end


function test_3d(n::Int; seed::Int)
    Random.seed!(seed)
    n_dims = 3
    test_points = initialize_vertex(n, n_dims=n_dims)
    tree = initialize_tree_3d(test_points)

    for point in test_points
        insert_point(tree, point, n_dims=n_dims)
    end

    # x,y,z = plot_simplex_3d(tree.simplices[1], tree.vertices)
    # p = plot3d(x, y, z, label="Points", size=(800, 800))
    # for i in 2:length(tree.simplices)
    #     if !tree.simplices[i].dead && all(tree.simplices[i].vertices.>0)
    #         x,y,z = plot_simplex_3d(tree.simplices[i], tree.vertices)
    #         plot3d!(x, y, z, label="Points", size=(800, 800))
    #     end
    # end

    # scatter3d!([x for x in map(x -> x.position[1], test_points)], [y for y in map(x -> x.position[2], test_points)], [z for z in map(x -> x.position[3], test_points)], label="Points", color=["red","blue","green", "yellow", "black"])

    # center, R = circumsphere(1, tree)
    # surface!(sphere(center, R), color=:red, alpha=0.3)
    check_delaunay(tree, n_dims=3)
    return tree
end

using Profile
using ProfileView
@ProfileView.profview test_3d(20,seed=1)