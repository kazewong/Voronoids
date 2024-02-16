function batch_locate!(output::AbstractArray, vertices::AbstractArray, tree::DelaunayTree)
    for i in 1:length(vertices)
        output[i] = locate(Vector{Int}(), vertices[i], tree)
    end
    output
end

function parallel_locate(vertices::Vector{Vector{Float64}}, tree::DelaunayTree)::Vector{Vector{Int}}
    output = Vector{Vector{Int}}(undef, length(vertices))
    batch_size = length(vertices) ÷ Threads.nthreads()
    Threads.@threads for i in 1:Threads.nthreads()
        if i == Threads.nthreads()
            slice = @view output[(i-1)*batch_size+1:end]
            batch_locate!(slice, vertices[(i-1)*batch_size+1:end], tree)
        else
            slice = @view output[(i-1)*batch_size+1:i*batch_size]
            batch_locate!(slice, vertices[(i-1)*batch_size+1:i*batch_size], tree)
        end
    end
    return output
end

function identify_nonconflict_points(vertices::Vector{Vector{Float64}}, tree::DelaunayTree)::Tuple{Vector{Vector{Int}},Vector{Int}}
    site_list = parallel_locate(vertices, tree)
    length_site_list = length(site_list)
    output = Vector{Bool}(undef, length_site_list)
    Threads.@threads for i in 1:length_site_list-1
        output[i] = all(isempty.(map(x->intersect(reduce(vcat,tree.neighbors_relation[site_list[i]]), reduce(vcat, tree.neighbors_relation[site_list[x]])), i+1:length_site_list)))
        # output[i] = all(isempty.(map(x -> intersect(site_list[i], site_list[x]), i+1:length_site_list)))
    end
    return site_list, output
end

# Try the sphere check for identifying nonconflict point later.

function make_new_neighbors(simplices::Vector{Vector{Int}}; n_dims::Int=3)::Vector{Tuple{Int, Int}}
    length_simplices = length(simplices)
    output = Vector{Tuple{Int, Int}}()
    for i in 1:length_simplices
        for j in i+1:length_simplices
            facet = common_facet(simplices[i], simplices[j])
            if length(facet) == n_dims
                push!(output, (i, j))
                push!(output, (j, i))
            end
        end
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

function batch_insert_point(vertices::Vector{Vector{Float64}}, tree::DelaunayTree)
    sites, conflict = identify_nonconflict_points(vertices, tree)
    sites = sites[conflict]
    vertices = vertices[conflict]
    vertex_ids = [length(tree.vertices) + i for i in 1:length(vertices)]
    simplices = Vector{Vector{Vector{Int}}}(undef, length(vertices))
    centers = Vector{Vector{Vector{Float64}}}(undef, length(vertices))
    radii = Vector{Vector{Float64}}(undef, length(vertices))
    neighbors_id = Vector{Vector{Tuple{Int, Int}}}(undef, length(vertices))
    new_neighbors_id = Vector{Vector{Tuple{Int, Int}}}(undef, length(vertices))
    vertices_simplex = Vector{Vector{Vector{Int}}}(undef, length(vertices))

    Threads.@threads for i in 1:length(vertices)
        _simplices, _centers, _radii, _neighbors_id, _new_neighbors_id = new_simplex(sites[i], vertices[i], ids[i], tree)
        simplices[i] = _simplices
        centers[i] = _centers
        radii[i] = _radii
        neighbors_id[i] = _neighbors_id
        new_neighbors_id[i] = _new_neighbors_id
    end
    return simplices, centers, radii, neighbors_id, new_neighbors_id
    # push!(tree.simplices, simplices...)
    # push!(tree.centers, centers...)
    # push!(tree.radii, radii...)
end

export parallel_locate, batch_locate!, identify_nonconflict_points, new_simplex, batch_insert_point