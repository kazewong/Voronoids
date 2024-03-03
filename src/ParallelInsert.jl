

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

function identify_conflicts(vertices::Vector{Vector{Float64}}, tree::DelaunayTree)::Tuple{Vector{Vector{Int}},Vector{Vector{Int}}}
    site_list = parallel_locate(vertices, tree)
    neighbor_list = Vector{Vector{Int}}(undef,length(site_list))
    Threads.@threads for i in 1:length(site_list)
        neighbor_list[i] = unique(mapreduce(x->tree.neighbors_relation[x], vcat, site_list[i]))
    end
    return site_list, neighbor_list
end

function add_to_occupancy!(occupancy::Dict{Int, Vector{Int}}, sites::Vector{Vector{Int}}, neighbors::Vector{Vector{Int}}, ids::Vector{Int})
    for i in 1:length(sites) # This might be able to be parallelized
        for neighbor in neighbors[i]
            if !haskey(occupancy, neighbor)
                occupancy[neighbor] = Vector{Int}()
            end
            occupancy[neighbor] = push!(occupancy[neighbor], ids[i])
        end
    end
end

function find_conflict_group(result::Vector{Int}, neighbors::Vector{Vector{Int}}, vertex_id::Int, occupancy::Dict{Int, Vector{Int}})::Vector{Int}
    push!(result, vertex_id)
    conflict_list = unique(mapreduce(x->occupancy[x], vcat, neighbors[vertex_id % (length(neighbors)+1)]))
    for conflict in conflict_list
        if conflict ∉ result
            result = find_conflict_group(result, neighbors, conflict, occupancy)
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
        conflict_group = find_conflict_group(Vector{Int}(), neighbors, vertex_id, occupancy)
        push!(output, conflict_group)
        push!(checked, conflict_group...)
    end
    return output
end

function queue_multiple_points!(channel::Channel{Tuple{Int, Vector{Float64},Vector{Int}}}, points::Vector{Vector{Float64}}, occupancy::Dict{Int, Vector{Int}},tree::DelaunayTree, lk:: ReentrantLock; batch_size::Int=256)
    partition = Iterators.partition(1:length(points), batch_size)
    println("Size of partition: ", length(points))
    for chunk in partition
        site_list, neighbor_list = identify_conflicts(points[chunk], tree)
        lock(lk)
        add_to_occupancy!(occupancy, site_list, neighbor_list, collect(chunk))
        unlock(lk)
        for i in 1:length(chunk)
            put!(channel, (chunk[i], points[chunk[i]], neighbor_list[i]))
        end
    end
    # put!(channel, (-1, [-1.], [-1]))
    println("Done queuing")
end

function consume_point!(id::Int, vertex::Vector{Float64}, neighbors::Vector{Int}, tree::DelaunayTree, lk::ReentrantLock, occupancy::Dict{Int, Vector{Int}}; n_dims::Int)
    lock(lk)
        add_vertex!(tree, vertex, n_dims=n_dims)
        for i in 1:length(neighbors)
            popfirst!(occupancy[neighbors[i]])
            if isempty(occupancy[neighbors[i]])
                delete!(occupancy, neighbors[i])
            end
        end
    unlock(lk)
end

function consume_multiple_points!(channel::Channel{Tuple{Int, Vector{Float64}, Vector{Int}}}, tree::DelaunayTree, queuing::Task, occupancy::Dict{Int, Vector{Int}},lk::ReentrantLock, n_dims::Int)
    wait_queue = Vector{Tuple{Int, Vector{Float64}, Vector{Int}}}()
    tasks = Vector{Tuple{Int, Task}}()
    while !istaskdone(queuing) || !isempty(channel) || !isempty(wait_queue)
        if !isempty(channel)
            points = take!(channel)
            # if points[1] == -1
            #     break
            # end
            if any(length.(map(x->occupancy[x],points[3])).!=1)
                push!(wait_queue, points)
            else
                push!(tasks,(points[1],Threads.@spawn consume_point!(points[1], points[2], points[3], tree, lk, occupancy, n_dims=n_dims)))
            end
        end
        if !isempty(wait_queue)
            for i in 1:length(wait_queue)
                if all(map(y->y[1]==wait_queue[i][1], map(x->occupancy[x],wait_queue[i][3])))
                    points = popat!(wait_queue, i)
                    push!(tasks,(points[1],Threads.@spawn consume_point!(points[1], points[2], points[3], tree, lk, occupancy, n_dims=n_dims)))
                end
            end
        end
    end
    println("Done consuming")
    println("Number of vertices: ", length(tree.vertices))
    println("is empty: ", isempty(channel))
    println("is empty: ", isempty(wait_queue))
    return tasks
end


function parallel_insert!(points::Vector{Vector{Float64}}, tree::DelaunayTree; n_dims::Int=3, batch_size::Int=256)
    update_channel = Channel{Tuple{Int, Vector{Float64}, Vector{Int}}}(batch_size+1)
    lk = ReentrantLock()
    occupancy = Dict{Int, Vector{Int}}()

    queue_multiple_points!(update_channel, points, occupancy, tree, lk)
    consume_points!(update_channel, tree, t1, occupancy, lk, n_dims)
    return update_channel, t1, t2
end

export parallel_locate, batch_locate, identify_conflicts, find_conflict_group, group_points, queue_multiple_points!, consume_multiple_points!, parallel_insert!