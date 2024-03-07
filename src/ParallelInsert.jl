

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

function queue_multiple_points!(channel::Channel{Tuple{Int, Vector{Float64},Vector{Int}}}, points::Vector{Vector{Float64}}, occupancy::Dict{Int, Vector{Int}},tree::DelaunayTree, lk:: ReentrantLock; batch_size::Int=256)
    partition = Iterators.partition(1:length(points), batch_size)
    println("Size of partition: ", length(points))
    for chunk in partition
        site_list, neighbor_list = identify_conflicts(points[chunk], tree)
        for i in 1:length(neighbor_list)
            neighbor_list[i] = unique(mapreduce(x->tree.neighbors_relation[x], vcat, neighbor_list[i]))
        end
        lock(lk)
        add_to_occupancy!(occupancy, site_list, neighbor_list, collect(chunk))
        unlock(lk)
        for i in 1:length(chunk)
            put!(channel, (chunk[i], points[chunk[i]], neighbor_list[i]))
        end
    end
    println("Done queuing")
end

function channel_to_queue(n_points::Int, channel::Channel{Tuple{Int, Vector{Float64},Vector{Int}}})::Vector{Tuple{Int, Vector{Float64},Vector{Int}}}
    queue = Vector{Tuple{Int, Vector{Float64},Vector{Int}}}()
    while !isempty(channel) || length(queue) < n_points
        push!(queue, take!(channel))
    end
    return queue
end

function find_placement(id::Int,neighbors::Vector{Int},  occupancy::Dict{Int, Vector{Int}})::Int
    return findfirst(x->x==id, unique(reduce(vcat, map(x->occupancy[x], neighbors))))
end

function find_placement(id::Vector{Int}, neighbors::Vector{Vector{Int}},  occupancy::Dict{Int, Vector{Int}})::Vector{Int}
    return map(x->find_placement(x[1], x[2], occupancy), zip(id, neighbors))
end

function add_multiple_vertex!(tree::DelaunayTree, vertices::Vector{Vector{Float64}}, lk::ReentrantLock; n_dims::Int)
    updates = Vector{TreeUpdate}(undef, length(vertices))
    n_vertex = length(tree.vertices)
    Threads.@threads for i in 1:length(vertices)
        updates[i] = make_update(n_vertex+i, vertices[i], tree, n_dims=n_dims)
    end
    for i in 1:length(vertices)
        insert_point!(tree, updates[i])
        add_point!(tree.kdtree, vertices[i])
    end
    return updates
end

function update_multiple_occupancy!(occupancy::Dict{Int, Vector{Int}}, neighbors::Vector{Vector{Int}}, ids::Vector{Int})
    Threads.@threads for neighbor in neighbors
        for i in 1:length(neighbor)
            popfirst!(occupancy[neighbor[i]])
            if isempty(occupancy[neighbor[i]])
                delete!(occupancy, neighbor[i])
            end
        end
    end
end

function consume_multiple_points!(wait_queue::Vector{Tuple{Int, Vector{Float64},Vector{Int}}}, tree::DelaunayTree, occupancy::Dict{Int, Vector{Int}}, n_dims::Int)
    timer = time()
    placement = find_placement(getindex.(wait_queue, 1), getindex.(wait_queue, 3), occupancy)

    println(maximum(placement))
    for i in 1:maximum(placement)
        index = findall(x->x==i, placement)
        println("Point inserted: $(length(index)), Point per second: $((length(index))/(time()-timer))")

        # println("Number of non-blocked points: ", sum(non_block_live_point))
        ids = getindex.(wait_queue[index], 1)
        vertices = getindex.(wait_queue[index], 2)
        neighbors = getindex.(wait_queue[index], 3)
        # add_multiple_vertex!(tree, vertices, lk, n_dims=n_dims)
        update_multiple_occupancy!(occupancy, neighbors, ids)
    end
end


function parallel_insert!(points::Vector{Vector{Float64}}, tree::DelaunayTree; n_dims::Int=3, batch_size::Int=256)
    update_channel = Channel{Tuple{Int, Vector{Float64}, Vector{Int}}}(batch_size+1)
    lk = ReentrantLock()
    occupancy = Dict{Int, Vector{Int}}()

    queue_multiple_points!(update_channel, points, occupancy, tree, lk)
    consume_points!(update_channel, tree, t1, occupancy, lk, n_dims)
    return update_channel, t1, t2
end

export parallel_locate, batch_locate, identify_conflicts, queue_multiple_points!, channel_to_queue, find_placement, consume_multiple_points!, parallel_insert!