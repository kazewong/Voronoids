

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

function identify_conflicts!(vertices::Vector{Vector{Float64}}, occupancy::Dict{Int, Vector{Int}}, tree::DelaunayTree)::Tuple{Vector{Vector{Int}},Vector{Vector{Int}}}
    #=
    As opposed to checking for conflict per pair of vertices, the idea is to use an occupancy list to indiciate whether a site has any conflict with its neighbors.

    This involves populating the occupancy list with the site indices, which has complexity O(n) as opposed to O(n^2) for the pair-wise check.

    The output contains the indices of vertices which occupies a certain simplex, which can be used to trace back to which vertices is a particular vertex conflicting with.

    This will benefit from pruning the dead simplices, which is not implemented yet.

    =#
    site_list = parallel_locate(vertices, tree)
    neighbor_list = Vector{Vector{Int}}(undef,length(site_list))
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
    return site_list, neighbor_list
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

function queue_multiple_points!(channel::Channel{Vector{Tuple{Vector{Float64}, Vector{Int}, Vector{Int}}}}, points::Vector{Vector{Float64}}, occupancy::Dict{Int, Vector{Int}},tree::DelaunayTree, lk:: ReentrantLock; batch_size::Int=256)
    partition = Iterators.partition(1:length(points), batch_size)
    for chunk in partition
        lock(lk)
        site_list, neighbor_list = identify_conflicts!(points[chunk], occupancy, tree)
        unlock(lk)
        groups = group_points(site_list, neighbor_list, occupancy)
        for i in 1:length(groups)
            output = Vector{Tuple{Vector{Float64}, Vector{Int}, Vector{Int}}}()
            for j in 1:length(groups[i])
                push!(output, (points[groups[i][j]], site_list[groups[i][j]], neighbor_list[groups[i][j]]))
            end
            put!(channel, output)
        end
    end
    put!(channel, [([-1.], [-1], [-1])])
    println("Done queuing")
end

function consume_points!(channel::Channel{Vector{Tuple{Vector{Float64}, Vector{Int}, Vector{Int}}}}, tree::DelaunayTree, queuing::Task, occupancy::Dict{Int, Vector{Int}},lock::ReentrantLock, n_dims::Int)
    while !istaskdone(queuing) || !isempty(channel)
        points = take!(channel)
        if points[1][2] == [-1.]
            break
        end
        if length(points) > 1
            serial_insert!(collect(map(x->x[1],points)), tree, n_dims=n_dims)
        else
            insert_point!(tree, make_update(points[1][1], points[1][2], tree, n_dims=n_dims))
        end
    end
end


function parallel_insert!(points::Vector{Vector{Float64}}, tree::DelaunayTree; n_dims::Int=3, batch_size::Int=256)
    update_channel = Channel{Vector{Tuple{Vector{Float64}, Vector{Int}, Vector{Int}}}}(batch_size+1)
    lk = ReentrantLock()
    occupancy = Dict{Int, Vector{Int}}()

    t1 = @async queue_multiple_points!(update_channel, points, occupancy, tree, lk)
    t2 = @async consume_points!(update_channel, tree, t1, occupancy, lk, n_dims)
    return update_channel, t1, t2
end

export parallel_locate, batch_locate, identify_conflicts!, find_conflict_group, group_points, queue_multiple_points!, consume_points!, parallel_insert!