using Revise
using BoundingSphere
using LinearAlgebra
# using Plots
# plotlyjs()
using Random
using TimerOutputs
include("Primitives.jl")
using Memoize

const tmr = TimerOutput()

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

@memoize function circumcircle(vectices::Vector{Vector{Float64}})
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

@memoize function circumsphere(vertices::Vector{Vector{Float64}}; n_dims::Int=3)
    if n_dims == 2
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
    elseif n_dims == 3
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
end

function initialize_tree_3d(positions::Vector{Vector{Float64}})::DelaunayTree
    center, radius = boundingsphere(positions)
    radius = radius*5
    first_vertex = center + [0, 0, radius]
    second_vertex = center + [radius * cos(0), radius * sin(0), -radius ]
    third_vertex = center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3), -radius]
    fourth_vertex = center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3), -radius]
    ghost_vertex = [center + [0, 0, radius], center + [0, 0, radius], center + [0, 0, radius], center + [radius  * cos(0), radius * sin(0), -radius]]
    verticies = [first_vertex, second_vertex, third_vertex, fourth_vertex, ghost_vertex...]

    id = [-7, -6, -5, -4, -3, -2, -1, 0]
    simplicies = [[1, 2, 3, 4], [5, 1, 2, 3], [6, 1, 3, 4], [7, 1, 4, 2], [8, 2, 3, 4]]
    dead = [false, false, false, false, false]
    centers = [center, [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]
    radii = [radius, 0, 0, 0, 0]

    parent_relation = [1, 2, 3, 4, 5]
    children_relation= Vector{Vector{Int}}([[],[],[],[],[]])
    step_children_relation = [Pair([],[]),Pair([],[]),Pair([],[]),Pair([],[]),Pair([],[])]
    neighbors_relation = [[2, 3, 4, 5], [1], [1], [1], [1]]
    return DelaunayTree(id, simplicies, dead, centers, radii, parent_relation, children_relation, step_children_relation, neighbors_relation, verticies)
end

function initialize_tree_2d(positions::Vector{Vector{Float64}})::DelaunayTree
    center, radius = boundingsphere(positions)
    radius = radius*3
    first_vertex = center + [radius * cos(0* pi / 3), radius * sin(0 * pi / 3)]
    second_vertex = center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3)]
    third_vertex = center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3)]
    ghost_vertex = [center + [radius * cos(0 * pi / 3), radius * sin(0 * pi / 3)], center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3)], center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3)]]
    verticies = [first_vertex, second_vertex, third_vertex, ghost_vertex...]

    id = [-5, -4, -3, -2, -1, 0]
    simplicies = [[1, 2, 3], [4, 1, 2], [6, 1, 3], [5, 2, 3]]
    dead = [false, false, false, false]
    centers = [center, [0, 0], [0, 0], [0, 0]]
    radii = [radius, 0, 0, 0]

    parent_relation = [1, 2, 3, 4]
    children_relation= Vector{Vector{Int}}([[],[],[],[]])
    step_children_relation = [Pair([],[]),Pair([],[]),Pair([],[]),Pair([],[])]
    neighbors_relation = [[2, 3, 4], [1], [1], [1]]
    return DelaunayTree(id, simplicies, dead, centers, radii, parent_relation, children_relation, step_children_relation, neighbors_relation, verticies)
end

function in_sphere(node_id::Int, point::Vector{Float64}, tree::DelaunayTree)::Bool
    return norm(point .- tree.centers[node_id]) < tree.radii[node_id]
end

function locate(visited_ids::Vector{Int}, output::Vector{Int}, vertex::Vector{Int}, current_node_id::Int, tree::DelaunayTree; n_dims::Int = 3)::Vector{Int}
    if current_node_id ∉ visited_ids && in_sphere(current_node_id, vertex, tree)
        if !tree.dead[current_node_id]
            push!(output, current_node_id)
            return output
        end
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

function locate(vertex::Vertex, tree::DelaunayTree; n_dims::Int = 3)::Vector{Int}
    alive_point_id = map(x->x.id, filter(x->!x.dead, collect(values(tree.simplices))))
    insphere = false
    node_id = 1
    while !insphere
        if in_sphere(alive_point_id[node_id], vertex, tree)
            insphere = true
            break
        else
            node_id += 1
        end
    end
    return find_all_neighbors(Vector{Int}(), alive_point_id[node_id], vertex, tree, n_dims=n_dims)
end

function find_all_neighbors(output::Vector{Int}, node_id::Int, point::Vertex, tree::DelaunayTree; n_dims=3)::Vector{Int}
    neighbors = tree.neighbors_relation[node_id]
    for neighbor_id in neighbors
        if neighbor_id ∉ output && in_sphere(neighbor_id, point, tree)
            push!(output, neighbor_id)
            find_all_neighbors(output, neighbor_id, point, tree)
        end
    end
    return output
end

function common_facet(simplex1::Vector{Int}, simplex2::Vector{Int}; n_dims::Int = 3)::Vector{Int}
    @timeit tmr "check intersect" common = intersect(simplex1, simplex2)
    if length(common) == n_dims
        return common
    else
        return []
    end
end

function insert_point(tree::DelaunayTree, point::Vertex; n_dims::Int=3)
    @timeit tmr "locating node" killed_nodes = locate(Vector{Int}(), Vector{Int}(), point, 1, tree, n_dims=n_dims)
    # println("killed_nodes: ", filter(x->tree.simplices[x].dead==false, killed_nodes))
    new_node_id = Vector{Int}()
    @timeit tmr "insert per killed nodes" for node_id in killed_nodes
        if !tree.dead[node_id].dead
            tree.dead[node_id].dead = true
            for neighbor_id in tree.neighbors_relation[node_id]
                # println(in_sphere(neighbor_id, point, tree, n_dims=n_dims))
                @timeit tmr "in sphere" if !in_sphere(neighbor_id, point, tree)
                    # println("neighbor_id: ", neighbor_id)
                    @timeit tmr "check facet" facet = common_facet(tree.vertices[tree.simplices[node_id]], tree.vertices[tree.simplices[neighbor_id]], n_dims=n_dims)
                    if length(facet) == n_dims
                        # Creating new node
                        new_id = tree.id[end] + 1
                        push!(tree.id, new_id)
                        push!(tree.simplices, [point.id, facet...])
                        push!(tree.dead, false)
                        center, radius = circumsphere([point.position, map(x->tree.vertices[x].position, facet)...], n_dims=n_dims)
                        push!(tree.centers, center)
                        push!(tree.radii, radius)

                        push!(tree.parent_relation, node_id)
                        push!(tree.children_relation, Vector{Int}())
                        push!(tree.step_children_relation, Vector{Pair(Vector{Int}, Vector{Int})}())
                        push!(tree.neighbors_relation, [neighbor_id])

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
                            tree.neighbors_relation[neighbor_id][killed_node_id] = new_node.id
                        end
                    end
                end
            end
        end
    end
    
    # println("len new_node_id: ", length(new_node_id))
    @timeit tmr "making new neighbor" for i in 1:length(new_node_id)
        for j in i+1:length(new_node_id)
            new_id1 = new_node_id[i]
            new_id2 = new_node_id[j]
            facet = common_facet(tree.vertices[tree.simplices[new_id1]], tree.vertices[tree.simplices[new_id2]], n_dims=n_dims)
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
                    if in_sphere(i, tree.vertices[j], tree) && j ∉ tree.simplices[i].vertices && all(tree.simplices[i].vertices.>0)
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
    @timeit tmr "initializing tree" tree = initialize_tree_2d(test_points)
    
    for point in test_points
        insert_point(tree, point, n_dims=n_dims)
    end
    
    # x,y = plot_simplex_2d(tree.simplices[1], tree.vertices)
    # plot(x, y, label="Points", size=(800, 800))
    # for i in 2:length(tree.simplices)
    #     if !tree.simplices[i].dead && all(tree.simplices[i].vertices.>0)
    #         x,y = plot_simplex_2d(tree.simplices[i], tree.vertices)
    #         plot!(x, y, label="Points", size=(800, 800))
    #     end
    # end
    
    # p = scatter!([x for x in map(x -> x.position[1], test_points)], [y for y in map(x -> x.position[2], test_points)], label="Points", color=["red","blue","green", "yellow", "black"])
    
    check_delaunay(tree, n_dims=2)

    return tree, p
end


function test_3d(n::Int; seed::Int)
    Random.seed!(seed)
    n_dims = 3
    @timeit tmr "initializing points" test_points = initialize_vertex(n, n_dims=n_dims)
    @timeit tmr "initializing tree" tree = initialize_tree_3d(test_points)

    for point in test_points
        @timeit tmr "insert_point" insert_point(tree, point, n_dims=n_dims)
    end

    # x,y,z = plot_simplex_3d(tree.simplices[1], tree.vertices)
    # p = plot3d(x, y, z, label="Points", size=(800, 800))
    # for i in 2:length(tree.simplices)
    #     if !tree.simplices[i].dead && all(tree.simplices[i].vertices.>0)
    #         x,y,z = plot_simplex_3d(tree.simplices[i], tree.vertices)
    #         plot3d!(x, y, z, label="Points", size=(800, 800))
    #     end
    # end

    # scatter3d!([x for x in map(x -> x.position[1], test_points)], [y for y in map(x -> x.position[2], test_points)], [z for z in map(x -> x.position[3], test_points)], label="Points", c=distinguishable_colors(n))

    # center, R = circumsphere(1, tree)
    # surface!(sphere(center, R), color=:red, alpha=0.3)
    # check_delaunay(tree, n_dims=3)
    return tree#, p
end


# @timeit tmr "test2d" tree, p = test_2d(10,seed=1234)
@timeit tmr "test3d" tree = test_3d(50,seed=1)
tmr
# display(p)