function batch_locate(vertices::AbstractArray, tree::DelaunayTree)
    output = Vector{Vector{Int}}(undef, length(vertices))
    for i in 1:length(vertices)
        output[i] = locate(Vector{Int}(), vertices[i], tree)
    end
    return output
end

function parallel_locate(vertices::Vector{Vector{Float64}}, tree::DelaunayTree)::Vector{Vector{Int}}
    output = Vector{Vector{Int}}(undef, length(vertices))
    chunks = collect(Iterators.partition(eachindex(vertices), max(length(vertices) ÷ Threads.nthreads(),1)))
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

function make_queue!(points::Vector{Vector{Float64}}, occupancy::Dict{Int, Vector{Int}},tree::DelaunayTree; batch_size::Int=256)::Vector{Tuple{Int, Vector{Float64},Vector{Int}}}
    queue = Vector{Tuple{Int, Vector{Float64},Vector{Int}}}(undef, length(points))
    partition = Iterators.partition(1:length(points), batch_size)
    println("Size of partition: ", length(points))
    for chunk in partition
        site_list, neighbor_list = identify_conflicts(points[chunk], tree)
        for i in 1:length(neighbor_list)
            neighbor_list[i] = unique(mapreduce(x->tree.neighbors_relation[x], vcat, neighbor_list[i]))
        end
        add_to_occupancy!(occupancy, site_list, neighbor_list, collect(chunk))
        for i in 1:length(chunk)
            queue[chunk[i]] = (chunk[i], points[chunk[i]], neighbor_list[i])
        end
    end
    println("Done queuing.")
    return queue
end

@enum EventState begin
    WAITING
    READY
    COMPLETED
    FAILED
end

mutable struct Event
    @atomic state::EventState
    id::Int
    blocked_by::Vector{Int}
    task::Task
    Event(tree::DelaunayTree, point::Vector{Float64}, id::Int, blocked_by::Vector{Int}; n_dims::Int = 3) = new(WAITING, id, blocked_by, Task(() -> add_vertex!(tree, point, n_dims=n_dims)))
end

function run_event(e::Event, event_list::Vector{Event},inserted::Vector{Bool})
    state = @atomic e.state
    if state == WAITING
        while state !== READY
            state, ok = @atomicreplace(e.state, state => READY)
            if ok
                if all(inserted[e.blocked_by]) || length(e.blocked_by) == 0
                    schedule(e.task)
                    break
                end    
            end
        end
        while state !== COMPLETED && state !== FAILED
            state, ok = @atomicreplace(e.state, state => COMPLETED)
            if ok
                inserted[e.id] = true
            end
            # if istaskdone(e.task)
            #     inserted[e.id] = true
            #     e.state = COMPLETED
                # for i in e.id+1:length(inserted)
                #     if event_list[i].state == WAITING
                #         run_event(event_list[i], event_list, inserted)
                #     end
                #     break
                # end
            # end
        end
    end
end

function make_event(id::Int, point::Vector{Float64}, neighbors::Vector{Int}, occupancy::Dict{Int, Vector{Int}}, tree::DelaunayTree; n_dims::Int=3)
    order_id = sort(unique(reduce(vcat, map(x->occupancy[x], neighbors))))
    blocked_by = order_id[1:findfirst(x->x==id, order_id)-1]
    return Event(tree, point, id, blocked_by, n_dims=n_dims)
end

function make_event(queue::Vector{Tuple{Int, Vector{Float64},Vector{Int}}}, occupancy::Dict{Int, Vector{Int}}, tree::DelaunayTree; n_dims::Int=3)
    events = Vector{Event}(undef, length(queue))
    Threads.@threads for i in 1:length(queue)
        events[i] = make_event(queue[i][1], queue[i][2], queue[i][3], occupancy, tree, n_dims=n_dims)
    end
    return events
end

function find_placement!(placement::Vector{Int}, start_id::Int, neighbors::Vector{Vector{Int}},  occupancy::Dict{Int, Vector{Int}})
    if placement[start_id] ==0
        order_id = sort(unique(reduce(vcat, map(x->occupancy[x], neighbors[start_id]))))
        last_id = findfirst(x->x==start_id, order_id)
        # println("Last id: ", last_id)
        if last_id == 1
            placement[start_id] = 1
        else
            for i in last_id-1:-1:1
                placement[start_id] = max(placement[start_id], find_placement!(placement, order_id[i], neighbors, occupancy) + 1)
            end
        end
    end
    return placement[start_id]
end

function find_placement(neighbors::Vector{Vector{Int}},  occupancy::Dict{Int, Vector{Int}})::Vector{Int}
    placement = fill(0, length(neighbors))
    start_id = length(neighbors)
    while !isnothing(findfirst(x->x==0, placement))
        start_id = findlast(x->x==0, placement)
        find_placement!(placement, start_id, neighbors, occupancy)
    end
    println("Unique placement: ", maximum(placement))
    return placement
end

function add_multiple_vertex!(tree::DelaunayTree, vertices::Vector{Vector{Float64}}; n_dims::Int)
    updates = Vector{TreeUpdate}(undef, length(vertices))
    n_vertex = length(tree.vertices)
    Threads.@threads for i in 1:length(vertices)
        updates[i] = make_update(i+n_vertex, vertices[i], tree, n_dims=n_dims)
    end
    for i in 1:length(vertices)
        insert_point!(tree, updates[i])
        add_point!(tree.kdtree, vertices[i])
    end
    return updates
end

function consume_multiple_points!(wait_queue::Vector{Tuple{Int, Vector{Float64},Vector{Int}}}, tree::DelaunayTree, occupancy::Dict{Int, Vector{Int}},n_dims::Int)
    timer = time()
    placement = find_placement(getindex.(wait_queue, 3), occupancy)
    all_vertices = getindex.(wait_queue, 2)
    println("Placement found: $(time()-timer)")

    timer = time()
    for i in 1:maximum(placement)
        index = findall(x->x==i, placement)
        # println("Number of non-blocked points: ", sum(non_block_live_point))
        vertices = all_vertices[index]
        add_multiple_vertex!(tree, vertices, n_dims=n_dims)
    end
    println("Point inserted: $(length(wait_queue)), Point per second: $((length(wait_queue))/(time()-timer))")

end


function parallel_insert!(points::Vector{Vector{Float64}}, tree::DelaunayTree; n_dims::Int=3, batch_size::Int=256)
    for chunk in Iterators.partition(1:length(points), batch_size)
        occupancy = Dict{Int, Vector{Int}}()
        queue = make_queue!(points[chunk], occupancy, tree)
        consume_multiple_points!(queue, tree, occupancy, n_dims)
    end

end

export parallel_locate, batch_locate, identify_conflicts, make_queue!, find_placement!, find_placement, consume_multiple_points!, parallel_insert!, Event, EventState, run_event, make_event