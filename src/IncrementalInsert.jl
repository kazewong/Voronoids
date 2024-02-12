using LinearAlgebra

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

function find_all_neighbors(output::Vector{Int}, node_id::Int, point::Vector{Float64}, tree::DelaunayTree; n_dims=3)::Vector{Int}
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

function locate(output::Vector{Int}, vertex::Vector{Float64}, tree::DelaunayTree; n_dims::Int = 3)::Vector{Int}
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
    common = intersect(simplex1, simplex2)
    if length(common) == n_dims
        return common
    else
        return []
    end
end

function insert_point(tree::DelaunayTree, point::Vector{Float64}; n_dims::Int=3)
    killed_nodes = locate(Vector{Int}(), point, tree, n_dims=n_dims)
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
    
    # println("len new_node_id: ", length(new_node_id))
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

export check_delaunay, in_sphere, find_all_neighbors, find_nearest_simplex, locate, common_facet, insert_point