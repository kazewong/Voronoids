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
        for i in keys(tree.simplices)
            for j in 1:length(tree.vertices)
                if in_sphere(i, tree.vertices[j], tree) && tree.vertices[j] ∉ tree.vertices[tree.simplices[i]] && all(tree.simplices[i].>8)
                    println("Error, point ", j, " is inside the circumcircle of simplex ", i)
                end
            end
        end
    elseif n_dims==2
        for i in keys(tree.simplices)
            for j in tree.vertices
                if in_sphere(i, j, tree) && j ∉ tree.vertices[tree.simplices[i]] && all(tree.simplices[i].>6)
                    println("Error, point ", j, " is inside the circumcircle of simplex ", i)
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
    return locate(output, simplex_id, vertex, tree)
end

function locate(output::Vector{Int}, simplex_id::Vector{Int}, vertex::Vector{Float64}, tree::DelaunayTree)::Vector{Int}
    for id in simplex_id
        if in_sphere(id, vertex, tree)
            push!(output, id)
            find_all_neighbors(output, id, vertex, tree)
        end
    end
    return unique(output)
end

function get_new_simplices(site::Int, vertex::Vector{Float64}, vertex_id::Int, tree::DelaunayTree; n_dims::Int=3)
    simplices = Vector{Vector{Int}}()
    centers = Vector{Vector{Float64}}()
    radii = Vector{Float64}()
    neighbors_id = Vector{Tuple{Int, Int}}()
    for neighbor_id in tree.neighbors_relation[site]
        if 270086 in tree.neighbors_relation[site]
            println("Neighbor id: ", neighbor_id)
            println("Site: ", site)
        end
        if !in_sphere(neighbor_id, vertex, tree)
            facet = common_facet(tree.simplices[site], tree.simplices[neighbor_id], n_dims=n_dims)
            if length(facet) == n_dims
                push!(simplices, [vertex_id, facet...])
                center, radius = circumsphere([vertex, tree.vertices[facet]...], n_dims=n_dims)
                push!(centers, center)
                push!(radii, radius)
                push!(neighbors_id, (neighbor_id, site))
            end
        end
    end
    return (simplices, centers, radii, neighbors_id)
end

function get_new_simplices(site::Int, vertex::Vector{Float64}, tree::DelaunayTree; n_dims::Int=3)
    vertex_id = length(tree.vertices) + 1
    return get_new_simplices(site, vertex, vertex_id, tree, n_dims=n_dims)
end


function make_new_neighbors(simplices::Vector{Vector{Int}}, simplices_id:: Vector{Int}; n_dims::Int=3)::Vector{Tuple{Int, Int}}
    length_simplices = length(simplices)
    output = Vector{Tuple{Int, Int}}()
    for i in 1:length_simplices
        for j in i+1:length_simplices
            facet = common_facet(simplices[i], simplices[j], n_dims=n_dims)
            if length(facet) == n_dims
                push!(output, (simplices_id[i], simplices_id[j]))
                push!(output, (simplices_id[j], simplices_id[i]))
            end
        end
    end
    return output
end


function make_update(id::Int, point::Vector{Float64}, killed_sites:: Vector{Int}, tree::DelaunayTree; n_dims::Int=3)::TreeUpdate
    if length(killed_sites) == 0
        println("Point: ", point)
        throw(ArgumentError("The point is already in the Delaunay triangulation"))
    end
    simplices = Vector{Vector{Vector{Int}}}(undef, length(killed_sites))
    simplices_ids = Vector{Vector{Int}}(undef, length(killed_sites))
    centers = Vector{Vector{Vector{Float64}}}(undef, length(killed_sites))
    radii = Vector{Vector{Float64}}(undef, length(killed_sites))
    neighbors_id = Vector{Vector{Tuple{Int, Int}}}(undef, length(killed_sites))
    simplices_counter = 1
    for i in 1:length(killed_sites)
        simplices[i], centers[i], radii[i], neighbors_id[i] = get_new_simplices(killed_sites[i], point, id, tree, n_dims=n_dims)
        if 270086 in map(x->x[1], neighbors_id[i])
            println("Vertex id: ", id)
            println("Site: ", killed_sites)
        end
        simplices_ids[i] = collect(simplices_counter:simplices_counter+length(simplices[i])-1)
        simplices_counter += length(simplices[i])
    end
    simplices = vcat(simplices...)
    simplices_ids = vcat(simplices_ids...)
    centers = vcat(centers...)
    radii = vcat(radii...)
    neighbors_id = vcat(neighbors_id...)
    new_neighbors_id = make_new_neighbors(simplices, simplices_ids, n_dims=n_dims)
    return TreeUpdate(point, killed_sites, simplices, simplices_ids, centers, radii, neighbors_id, new_neighbors_id)
end

function make_update(id::Int, point::Vector{Float64}, tree::DelaunayTree; n_dims::Int=3)::TreeUpdate
    killed_sites = locate(Vector{Int}(), point, tree)
    return make_update(id, point, killed_sites, tree, n_dims=n_dims)
end

function make_update(point::Vector{Float64}, tree::DelaunayTree; n_dims::Int=3)::TreeUpdate
    killed_sites = locate(Vector{Int}(), point, tree)
    id = length(tree.vertices) + 1
    return make_update(id, point, killed_sites, tree, n_dims=n_dims)
end

function add_vertex!(tree::DelaunayTree, point::Vector{Float64}; n_dims::Int=3)
    update = make_update(point, tree, n_dims=n_dims)
    insert_point!(tree, update)
end

function serial_insert!(points::Vector{Vector{Float64}}, tree::DelaunayTree; n_dims::Int=3)
    for point in points
        add_vertex!(tree, point, n_dims=n_dims)
    end
end

export circumsphere, check_delaunay, in_sphere, find_all_neighbors, find_nearest_simplex, locate, common_facet, get_new_simplices, make_new_neighbors, make_update, add_vertex!