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

function plot_simplex_2d(simplex_Id::Int, tree::DelaunayTree)
    x, y = [], []
    vertices = tree.vertices[tree.simplices[simplex_Id]]
    for vertex_id in [1,2,3,1]
        push!(x, vertices[vertex_id][1])
        push!(y, vertices[vertex_id][2])
    end
    return x, y
end

function plot_simplex_3d(simplex_Id::Int, tree::KDDelaunayTree)
    x, y, z = [], [], []
    vertices = tree.vertices[tree.simplices[simplex_Id]]
    for vertex_id in [1,2,3,1,4,3,2,4]
        push!(x, vertices[vertex_id][1])
        push!(y, vertices[vertex_id][2])
        push!(z, vertices[vertex_id][3])
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

@timeit tmr "circumsphere" function circumsphere(vertices::Vector{Vector{Float64}}; n_dims::Int=3)
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

        return ([X, Y], R)
    elseif n_dims == 3
        v1 = vertices[1]
        v2 = vertices[2]
        v3 = vertices[3]
        v4 = vertices[4]

        if (v1==v2) || (v1==v3) || (v1==v4) || (v2==v3) || (v2==v4) || (v3==v4)
            return ((0, 0, 0), 0)
        end

        @timeit tmr "make column" length_column = [
            v1[1]^2 + v1[2]^2 + v1[3]^2;
            v2[1]^2 + v2[2]^2 + v2[3]^2;
            v3[1]^2 + v3[2]^2 + v3[3]^2;
            v4[1]^2 + v4[2]^2 + v4[3]^2;
        ]

        @timeit tmr "make a" a = det([v1[1] v1[2] v1[3] 1;
            v2[1] v2[2] v2[3] 1;
            v3[1] v3[2] v3[3] 1;
            v4[1] v4[2] v4[3] 1])

        @timeit tmr "make Dx" Dx = det([length_column[1] v1[2] v1[3] 1;
            length_column[2] v2[2] v2[3] 1;
            length_column[3] v3[2] v3[3] 1;
            length_column[4] v4[2] v4[3] 1])

        @timeit tmr "make Dy" Dy = - det([length_column[1] v1[1] v1[3] 1;
            length_column[2] v2[1] v2[3] 1;
            length_column[3] v3[1] v3[3] 1;
            length_column[4] v4[1] v4[3] 1])

        @timeit tmr "make Dz" Dz = det([length_column[1] v1[1] v1[2] 1;
            length_column[2] v2[1] v2[2] 1;
            length_column[3] v3[1] v3[2] 1;
            length_column[4] v4[1] v4[2] 1])

        @timeit tmr "make center" center = [Dx/2/a, Dy/2/a, Dz/2/a]
        @timeit tmr "make radius" radius = sqrt((v1[1]-center[1])^2 + (v1[2]-center[2])^2 + (v1[3]-center[3])^2)

        return ([Dx/2/a,Dy/2/a,Dz/2/a], radius) # Return the center coordinates and the radius
    end
end

function initialize_tree_3d(positions::Vector{Vector{Float64}})::KDDelaunayTree
    center, radius = boundingsphere(positions)
    radius = radius*5
    first_vertex = center + [0, 0, radius]
    second_vertex = center + [radius * cos(0), radius * sin(0), -radius ]
    third_vertex = center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3), -radius]
    fourth_vertex = center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3), -radius]
    ghost_vertex = [center + [0, 0, radius], center + [0, 0, radius], center + [0, 0, radius], center + [radius  * cos(0), radius * sin(0), -radius]]
    vertices = [first_vertex, second_vertex, third_vertex, fourth_vertex, ghost_vertex...]

    simplicies = [[1, 2, 3, 4], [5, 1, 2, 3], [6, 1, 3, 4], [7, 1, 4, 2], [8, 2, 3, 4]]
    vertices_simplex = [[1,2,3,4], [1, 2, 4, 5], [1, 2, 3, 5], [1, 3, 4, 5], [2], [3], [4], [5]]
    centers = [center, [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]
    radii = [radius, 0, 0, 0, 0]
    dead = [false, false, false, false, false]


    neighbors_relation = [[2, 3, 4, 5], [1], [1], [1], [1]]
    kd_tree = KDTree(reduce(hcat, vertices))
    return KDDelaunayTree(vertices, dead, kd_tree, simplicies, vertices_simplex, centers, radii, neighbors_relation)
end

function initialize_tree_2d(positions::Vector{Vector{Float64}})::DelaunayTree
    center, radius = boundingsphere(positions)
    radius = radius*3
    first_vertex = center + [radius * cos(0* pi / 3), radius * sin(0 * pi / 3)]
    second_vertex = center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3)]
    third_vertex = center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3)]
    ghost_vertex = [center + [radius * cos(0 * pi / 3), radius * sin(0 * pi / 3)], center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3)], center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3)]]
    vertices = [first_vertex, second_vertex, third_vertex, ghost_vertex...]

    simplicies = [[1, 2, 3], [4, 1, 2], [6, 1, 3], [5, 2, 3]]
    centers = [center, [0, 0], [0, 0], [0, 0]]
    radii = [radius, 0, 0, 0]
    dead = [false, false, false, false]

    neighbors_relation = [[2, 3, 4], [1], [1], [1]]
    kd_tree = KDTree(reduce(hcat, vertices))
    return KDDelaynayTree(vertices, dead, kd_tree, simplicies, centers, radii, neighbors_relation)
end

function in_sphere(node_id::Int, point::Vector{Float64}, tree::KDDelaunayTree)::Bool
    return norm(point .- tree.centers[node_id]) < tree.radii[node_id]
end

function find_all_neighbors(output::Vector{Int}, node_id::Int, point::Vector{Float64}, tree::KDDelaunayTree; n_dims=3)::Vector{Int}
    neighbors = tree.neighbors_relation[node_id]
    for neighbor_id in neighbors
        if neighbor_id ∉ output && in_sphere(neighbor_id, point, tree)
            push!(output, neighbor_id)
            find_all_neighbors(output, neighbor_id, point, tree)
        end
    end
    return output
end

function find_nearest_simplex(point::Vector{Float64}, tree::KDDelaunayTree)::Vector{Int}
    return tree.vertices_simplex[nn(tree.kdtree, point)[1]]
end

function locate(output::Vector{Int}, vertex::Vector{Float64}, tree::KDDelaunayTree; n_dims::Int = 3)::Vector{Int}
    simplex_id = find_nearest_simplex(vertex, tree)
    simplex_id = filter(x->tree.dead[x]==false, simplex_id)
    for id in simplex_id
        if in_sphere(id, vertex, tree)
            push!(output, id)
            find_all_neighbors(output, id, vertex, tree)
        end
    end
    return unique(output)
end


function common_facet(simplex1::Vector{Int}, simplex2::Vector{Int}; n_dims::Int = 3)::Vector{Int}
    @timeit tmr "check intersect" common = intersect(simplex1, simplex2)
    if length(common) == n_dims
        return common
    else
        return []
    end
end

function insert_point(tree::KDDelaunayTree, point::Vector{Float64}; n_dims::Int=3)
    @timeit tmr "locating node" killed_nodes = locate(Vector{Int}(), point, tree, n_dims=n_dims)
    new_node_id = Vector{Int}()
    push!(tree.vertices_simplex, Vector{Int}())
    vertex_id = length(tree.vertices) + 1
    @timeit tmr "insert per killed nodes" for node_id in killed_nodes
        if !tree.dead[node_id]
            tree.dead[node_id] = true
            for neighbor_id in tree.neighbors_relation[node_id]
                if !in_sphere(neighbor_id, point, tree)
                    facet = common_facet(tree.simplices[node_id], tree.simplices[neighbor_id], n_dims=n_dims)
                    if length(facet) == n_dims
                        # Creating new node
                        new_id = length(tree.simplices) + 1
                        push!(new_node_id, new_id)
                        push!(tree.simplices, [vertex_id, facet...])
                        push!(tree.dead, false)
                        @timeit tmr "circumsphere" center, radius = circumsphere([point, tree.vertices[facet]...], n_dims=n_dims)
                        push!(tree.centers, center)
                        push!(tree.radii, radius)
                        push!(tree.neighbors_relation, [neighbor_id])

                        push!(tree.vertices_simplex[vertex_id], new_id)
                        for i in facet
                            push!(tree.vertices_simplex[i], new_id)
                        end

                        # Updating neighbor relationship for the neighbor of the killed node
                        killed_node_id = findfirst(x->x==node_id, tree.neighbors_relation[neighbor_id])
                        if killed_node_id !== nothing
                            tree.neighbors_relation[neighbor_id][killed_node_id] = new_id
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
            facet = common_facet(tree.simplices[new_id1], tree.simplices[new_id2], n_dims=n_dims)
            if length(facet) == n_dims
                push!(tree.neighbors_relation[new_id1], new_id2)
                push!(tree.neighbors_relation[new_id2], new_id1)
            end
        end
    end
    push!(tree.vertices,point)
    add_point!(tree.kdtree, point)
end

function check_delaunay(tree::KDDelaunayTree; n_dims::Int=3)
    if n_dims==3
        for i in 1:length(tree.simplices)
            if !tree.dead[i]
                for j in 1:length(tree.vertices)
                    if in_sphere(i, tree.vertices[j], tree) && tree.vertices[j] ∉ tree.vertices[tree.simplices[i]] && all(tree.simplices[i].>8)
                        println("Error, point ", j, " is inside the circumcircle of simplex ", i)
                    end
                end
            end
        end
    elseif n_dims==2
        for i in 1:length(tree.simplices)
            if !tree.dead[i]
                for j in tree.vertices
                    if in_sphere(i, j, tree) && j ∉ tree.vertices[tree.simplices[i]] && all(tree.simplices[i].>6)
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
    test_points = [rand(n_dims) for i in 1:n]
    @timeit tmr "initializing tree" tree = initialize_tree_2d(test_points)
    
    for point in test_points
        insert_point(tree, point, n_dims=n_dims)
    end
    
    x,y = plot_simplex_2d(1, tree)
    plot(x, y, label="Points", size=(800, 800))
    for i in 2:length(tree.simplices)
        if !tree.dead[i] && all(tree.simplices[i].>6)
            x,y = plot_simplex_2d(i, tree)
            plot!(x, y, label="Points", size=(800, 800))
        end
    end
    
    p = scatter!(map(x -> x[1], test_points), map(x -> x[2], test_points), label="Points", c=distinguishable_colors(n)) 
    check_delaunay(tree, n_dims=2)

    return tree, p
end


function test_3d(n::Int; seed::Int)
    Random.seed!(seed)
    n_dims = 3
    @timeit tmr "initializing points" test_points = [rand(n_dims) for i in 1:n]
    @timeit tmr "initializing tree" tree = initialize_tree_3d(test_points)

    for point in test_points
        @timeit tmr "insert_point" insert_point(tree, point, n_dims=n_dims)
    end

    # x,y,z = plot_simplex_3d(1, tree)
    # plot3d(x, y, z, label="Points", size=(800, 800))
    # for i in 2:length(tree.simplices)
    #     if !tree.dead[i] && all(tree.simplices[i].>8)
    #         x,y,z = plot_simplex_3d(i, tree)
    #         plot3d!(x, y, z, label="Points", size=(800, 800))
    #     end
    # end

    # p = scatter3d!(map(x -> x[1], test_points), map(x -> x[2], test_points), map(x -> x[3], test_points), label="Points", c=distinguishable_colors(n))
    # check_delaunay(tree, n_dims=3)
    return tree
end


# @timeit tmr "test2d" tree, p = test_2d(2,seed=1234)
# @timeit tmr "test3d" tree = test_3d(500000,seed=1)
# tmr

function parallel_locate(tree::KDDelaunayTree, points::Vector{Vector{Float64}}, tree::KDDelaunayTree; n_dims::Int=3)::Vector{Vector{Int}}
    Threads.@threads for point in points
        locate(Vector{Int}(), point, tree, n_dims=n_dims)
    end
end

function parallelInsert(tree::DelaunayTree, points::Vector{Vector{Float64}})
    # for point in points
    #     insert_point(tree, point, n_dims=3)
    # end
end