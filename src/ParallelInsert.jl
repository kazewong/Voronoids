

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

function identify_conflicts(vertices::Vector{Vector{Float64}}, tree::DelaunayTree)::Tuple{Vector{Vector{Int}},Vector{Vector{Int}},Dict{Int,Vector{Int}}}
    #=
    As opposed to checking for conflict per pair of vertices, the idea is to use an occupancy list to indiciate whether a site has any conflict with its neighbors.

    This involves populating the occupancy list with the site indices, which has complexity O(n) as opposed to O(n^2) for the pair-wise check.

    The output contains the indices of vertices which occupies a certain simplex, which can be used to trace back to which vertices is a particular vertex conflicting with.

    This will benefit from pruning the dead simplices, which is not implemented yet.

    =#
    site_list = parallel_locate(vertices, tree)
    neighbor_list = Vector{Vector{Int}}(undef,length(site_list))
    occupancy = Dict{Int,Vector{Int}}()
    Threads.@threads for i in 1:length(site_list)
        neighbor_list[i] = unique(mapreduce(x->tree.neighbors_relation[x], vcat, site_list[i]))
    end
    for i in 1:length(site_list) # This might be able to be parallelized
        for neighbor in neighbor_list[i]
            if !haskey(occupancy, neighbor)
                occupancy[neighbor] = Vector{Int}()
            end
            occupancy[neighbor] = push!(occupancy[neighbor], i)
        end
    end
    return site_list, neighbor_list, occupancy
end

function find_conflict_group(result::Vector{Int}, neighbors::Vector{Vector{Int}}, occupancy::Dict{Int, Vector{Int}},vertex_id::Int)::Vector{Int}
    push!(result, vertex_id)
    conflict_list = unique(mapreduce(x->occupancy[x], vcat, neighbors[vertex_id]))
    for conflict in conflict_list
        if conflict ∉ result
            result = find_conflict_group(result, neighbors, occupancy, conflict)
        end
    end
    return result
end

function group_points(site_list::Vector{Vector{Int}}, neighbors::Vector{Vector{Int}}, occupancy::Dict{Int, Vector{Int}})::Vector{Vector{Int}}
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


export parallel_locate, batch_locate, identify_conflicts, find_conflict_group, group_points