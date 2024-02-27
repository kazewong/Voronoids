using LinearAlgebra

function circumsphere(vertices::Vector{Vector{Float64}}; n_dims::Int=3)
    #=
    This version of the code needs to comptue 4 determinants, which doesn't seem to be ideal to me.
    There should be an easier way to determine the centers and the radius.
    =#
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

        length_column = [
            v1[1]^2 + v1[2]^2 + v1[3]^2;
            v2[1]^2 + v2[2]^2 + v2[3]^2;
            v3[1]^2 + v3[2]^2 + v3[3]^2;
            v4[1]^2 + v4[2]^2 + v4[3]^2;
        ]

        a = det([v1[1] v1[2] v1[3] 1;
            v2[1] v2[2] v2[3] 1;
            v3[1] v3[2] v3[3] 1;
            v4[1] v4[2] v4[3] 1])

        Dx = det([length_column[1] v1[2] v1[3] 1;
            length_column[2] v2[2] v2[3] 1;
            length_column[3] v3[2] v3[3] 1;
            length_column[4] v4[2] v4[3] 1])

        Dy = - det([length_column[1] v1[1] v1[3] 1;
            length_column[2] v2[1] v2[3] 1;
            length_column[3] v3[1] v3[3] 1;
            length_column[4] v4[1] v4[3] 1])

        Dz = det([length_column[1] v1[1] v1[2] 1;
            length_column[2] v2[1] v2[2] 1;
            length_column[3] v3[1] v3[2] 1;
            length_column[4] v4[1] v4[2] 1])

        center = [Dx/2/a, Dy/2/a, Dz/2/a]
        radius = sqrt((v1[1]-center[1])^2 + (v1[2]-center[2])^2 + (v1[3]-center[3])^2)

        return ([Dx/2/a,Dy/2/a,Dz/2/a], radius) # Return the center coordinates and the radius
    end
end

function check_delaunay(tree::DelaunayTree; n_dims::Int=3)
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

function in_sphere(node_id::Int, point::Vector{Float64}, tree::DelaunayTree)::Bool
    return norm(point .- tree.centers[node_id]) < tree.radii[node_id]
end

function find_all_neighbors(output::Vector{Int}, node_id::Int, point::Vector{Float64}, tree::DelaunayTree)::Vector{Int}
    neighbors = tree.neighbors_relation[node_id]
    for neighbor_id in neighbors
        if neighbor_id ∉ output && in_sphere(neighbor_id, point, tree)
            push!(output, neighbor_id)
            find_all_neighbors(output, neighbor_id, point, tree)
        end
    end
    return output
end

function find_nearest_simplex(point::Vector{Float64}, tree::DelaunayTree)::Vector{Int}
    return tree.vertices_simplex[nn(tree.kdtree, point)[1]]
end

function common_facet(simplex1::Vector{Int}, simplex2::Vector{Int}; n_dims::Int = 3)::Vector{Int}
    common = intersect(simplex1, simplex2)
    if length(common) == n_dims
        return common
    else
        return []
    end
end

function locate(output::Vector{Int}, vertex::Vector{Float64}, tree::DelaunayTree)::Vector{Int}
    simplex_id = find_nearest_simplex(vertex, tree)
    for id in simplex_id[tree.dead[simplex_id] .== false]
        if in_sphere(id, vertex, tree)
            push!(output, id)
            find_all_neighbors(output, id, vertex, tree)
        end
    end
    return unique(output)
end

function insert_point(tree::DelaunayTree, point::Vector{Float64}; n_dims::Int=3)
    killed_nodes = locate(Vector{Int}(), point, tree)
    new_node_id = Vector{Int}()
    push!(tree.vertices_simplex, Vector{Int}())
    vertex_id = length(tree.vertices) + 1
    for node_id in killed_nodes
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
                        center, radius = circumsphere([point, tree.vertices[facet]...], n_dims=n_dims)
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
    
    #=
    This part is O(n^2) and can be optimized.
    I think by checking the vertices and make use of the vertice to simplex assignment, this could be reduced to O(n)
    =#
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
    push!(tree.vertices,point)
    add_point!(tree.kdtree, point)
end

export check_delaunay, in_sphere, find_all_neighbors, find_nearest_simplex, locate, common_facet, insert_point, circumsphere