

function batch_locate(vertices::AbstractArray, tree::DelaunayTree)
    output = Vector{Vector{Int}}(undef, length(vertices))
    for i in 1:length(vertices)
        output[i] = locate(Vector{Int}(), vertices[i], tree)
    end
    return output
end

function parallel_locate(vertices::Vector{Vector{Float64}}, tree::DelaunayTree)::Vector{Vector{Int}}
    output = Vector{Vector{Int}}(undef, length(vertices))
    chunks = collect(Iterators.partition(eachindex(vertices), length(vertices) ÷ Threads.nthreads()))
    Threads.@threads for chunk in chunks
        output[chunk] = batch_locate(vertices[chunk], tree)
    end
    return output
end

function identify_conflicts(vertices::Vector{Vector{Float64}}, tree::DelaunayTree)::Tuple{Vector{Vector{Int}},Vector{Vector{Int}},Vector{Vector{Int}}}
    #=
    As opposed to checking for conflict per pair of vertices, the idea is to use an occupancy list to indiciate whether a site has any conflict with its neighbors.

    This involves populating the occupancy list with the site indices, which has complexity O(n) as opposed to O(n^2) for the pair-wise check.

    The output contains the indices of vertices which occupies a certain simplex, which can be used to trace back to which vertices is a particular vertex conflicting with.

    This will benefit from pruning the dead simplices, which is not implemented yet.

    =#
    site_list = parallel_locate(vertices, tree)
    neighbor_list = Vector{Vector{Int}}(undef,length(site_list))
    occupancy = Vector{Vector{Int}}(undef, length(tree.simplices))
    Threads.@threads for i in 1:length(site_list)
        neighbor_list[i] = filter(x->x ∉ site_list[i], unique(reduce(vcat, tree.neighbors_relation[site_list[i]])))
    end
    for i in 1:length(site_list) # This might be able to be parallelized
        for neighbor in neighbor_list[i]
            if !isdefined(occupancy,neighbor) 
                occupancy[neighbor] = Vector{Int}()
            end
            occupancy[neighbor] = push!(occupancy[neighbor], i)
        end
    end
    return site_list, neighbor_list, occupancy
end

function find_conflict_group(result::Vector{Int}, neighbors::Vector{Vector{Int}}, occupancy::Vector{Vector{Int}},vertex_id::Int)::Vector{Int}
    push!(result, vertex_id)
    conflict_list = unique(reduce(vcat, occupancy[neighbors[vertex_id]]))
    for conflict in conflict_list
        if conflict ∉ result
            result = find_conflict_group(result, neighbors, occupancy, conflict)
        end
    end
    return result
end

function group_points(site_list::Vector{Vector{Int}}, neighbors::Vector{Vector{Int}}, occupancy::Vector{Vector{Int}})::Vector{Vector{Int}}
    #=
    This function groups the vertices which occupies the same simplices together.
    =#
    output = Vector{Vector{Int}}()
    checked = Vector{Int}()
    while length(checked) < length(site_list)
        vertex_id = findfirst(x->x ∉ checked, 1:length(site_list))
        conflict_group = find_conflict_group(Vector{Int}(), neighbors, occupancy, vertex_id)
        push!(output, conflict_group)
        push!(checked, conflict_group...)
    end
    return output
end


function new_simplex(sites::Vector{Int}, vertex::Vector{Float64}, vertex_id::Int, tree::DelaunayTree; n_dims::Int=3)
    simplices = Vector{Vector{Int}}()
    centers = Vector{Vector{Float64}}()
    radii = Vector{Float64}()
    neighbors_id = Vector{Tuple{Int, Int}}()
    for site in sites
        if !tree.dead[site]
            neighbor_sites = tree.neighbors_relation[site]
            for neighbor_site in neighbor_sites
                if neighbor_site ∉ sites
                    facet = common_facet(tree.simplices[site], tree.simplices[neighbor_site], n_dims=n_dims)
                    push!(simplices, [vertex_id, facet...])
                    center, radius = circumsphere([vertex, tree.vertices[facet]...], n_dims=n_dims)
                    push!(centers, center)
                    push!(radii, radius)
                    push!(neighbors_id, (neighbor_site, site))
                end
            end
        end
    end
    new_neighbors_id = make_new_neighbors(simplices)
        
    return simplices, centers, radii, neighbors_id, new_neighbors_id
end

function batch_insert_point(vertices::Vector{Vector{Float64}}, tree::DelaunayTree; n_dims::Int=3)#::Vector{Bool}
    sites, conflict = identify_conflicts(vertices, tree)
    sites = sites[conflict]
    vertices = vertices[conflict]
    vertex_ids = [length(tree.vertices) + i for i in 1:length(vertices)]
    simplices = Vector{Vector{Vector{Int}}}(undef, length(vertices))
    centers = Vector{Vector{Vector{Float64}}}(undef, length(vertices))
    radii = Vector{Vector{Float64}}(undef, length(vertices))
    neighbors_id = Vector{Vector{Tuple{Int, Int}}}(undef, length(vertices))
    new_neighbors_id = Vector{Vector{Tuple{Int, Int}}}(undef, length(vertices))

    Threads.@threads for i in 1:length(vertices)
        _simplices, _centers, _radii, _neighbors_id, _new_neighbors_id = new_simplex(sites[i], vertices[i], vertex_ids[i], tree)
        simplices[i] = _simplices
        centers[i] = _centers
        radii[i] = _radii
        neighbors_id[i] = _neighbors_id
        new_neighbors_id[i] = _new_neighbors_id
    end
    simplices_id = [collect(length(tree.simplices)+1:length(tree.simplices)+length(simplices[1]))]
    for i in 2:length(simplices)
        push!(simplices_id, collect(simplices_id[i-1][end]+1:simplices_id[i-1][end]+length(simplices[i])))
    end

    new_neighbors_id = map(y->map(x->(simplices_id[y][x[1]],simplices_id[y][x[2]]), new_neighbors_id[y]),1:length(new_neighbors_id))

    simplices = reduce(vcat, simplices)
    centers = reduce(vcat, centers)
    radii = reduce(vcat, radii)
    neighbors_id = reduce(vcat, neighbors_id)
    new_neighbors_id = reduce(vcat, new_neighbors_id)
    simplices_id = reduce(vcat,simplices_id)#[length(tree.simplices) + i for i in 1:length(simplices)]
    simplices_id_repeated = reduce(vcat, fill.(simplices_id, n_dims+1))

    push!(tree.simplices, simplices...)
    push!(tree.dead, fill(false, length(simplices))...)
    push!(tree.centers, centers...)
    push!(tree.radii, radii...)

    Threads.@threads for i in 1:length(reduce(vcat, sites))
        tree.dead[reduce(vcat, sites)[i]] = true
    end

    # Replacing neighbor relationship
    push!(tree.neighbors_relation, collect(map(x->[x[1]], neighbors_id))...)
    for i in 1:length(neighbors_id)
        tree.neighbors_relation[neighbors_id[i][1]][tree.neighbors_relation[neighbors_id[i][1]].== neighbors_id[i][2]] .= simplices_id[i]
    end

    for i in 1:length(new_neighbors_id)
        push!(tree.neighbors_relation[new_neighbors_id[i][1]], new_neighbors_id[i][2])
    end

    # Updating vertex simplex relationship
    all_vertices = reduce(vcat, simplices)
    for i in 1:length(vertices)
        push!(tree.vertices_simplex, Vector{Int}())
    end
    for i in 1:length(all_vertices)
        tree.vertices_simplex[all_vertices[i]] = push!(tree.vertices_simplex[all_vertices[i]], simplices_id_repeated[i])
    end

    # Adding vertices
    push!(tree.vertices, vertices...)
    for vertex in vertices
        add_point!(tree.kdtree, vertex)
    end
    return conflict
end

function insert_point_parallel!(tree::DelaunayTree, point::Vector{Vector{Float64}}; n_dims::Int=3, batch_factor::Int=4)
    batch_size = Threads.nthreads() * batch_factor
    inserted = fill(false, length(point))
    while any(inserted .== false)
        current_points = findall(x->x==false, inserted)
        if length(current_points) > batch_size
            current_points = current_points[1:batch_size]
            conflict = batch_insert_point(point[current_points], tree, n_dims=n_dims)
            inserted[current_points[conflict]] .= true    
        else
            for point_id in current_points
                insert_point(tree, point[point_id], n_dims=n_dims)
                inserted[point_id] = true
            end
        end
    end
end

export parallel_locate, batch_locate, identify_conflicts, find_conflict_group, group_points, new_simplex, batch_insert_point, insert_point_parallel!