using Revise
using Voronoids
using TimerOutputs
using BenchmarkTools
using ThreadPinning
using Random

pinthreads(:cores)

const tmr = TimerOutput()

Random.seed!(1234)

n = 10000
n_dims = 3
batch_size = 10000

test_points = [rand(n_dims) for i in 1:n]
test_points2 = [rand(n_dims) for i in 1:1e5]
if n_dims == 3
    tree = initialize_tree_3d(test_points)
else
    tree = initialize_tree_2d(test_points)
end

for i in 1:n
    add_vertex!(tree, test_points[i], n_dims=n_dims)
end

parallel_insert!(test_points2, tree, n_dims=n_dims, batch_size=batch_size)

parallel_tree = deepcopy(tree)

test_points3 = [rand(n_dims) for i in 1:1e5]
inserted = fill(false, length(test_points3))

n_parallel = 10000

occupancy = Dict{Int, Vector{Int}}()
queue = make_queue!(test_points3, occupancy, parallel_tree)
placement = find_placement(getindex.(queue, 3), occupancy)
all_vertices = getindex.(queue, 2)
index = findall(x->x==2, placement)
vertices = all_vertices[index]
updates = Vector{TreeUpdate}(undef, length(vertices))

function batch_update(vertices, parallel_tree, n_dims)
    updates = Vector{TreeUpdate}(undef, length(vertices))
    for i in 1:length(vertices)
        updates[i] = make_update(i+length(parallel_tree.vertices), vertices[i], parallel_tree, n_dims=n_dims)
    end
    return updates
end

chunks = collect(Iterators.partition(1:length(vertices), max(length(vertices) ÷ Threads.nthreads(),1)))
Threads.@threads for chunk in chunks
    updates[chunk] = batch_update(vertices[chunk], parallel_tree, n_dims)
end

tasks = fetch.(map(chunks) do chunk
    Threads.@spawn updates[chunk] = batch_update(vertices[chunk], parallel_tree, n_dims)
end)

# parallel_insert!(test_points3, parallel_tree, n_dims=n_dims, batch_size=batch_size)

# channel = Channel{Tuple{Int,TreeUpdate}}(n_parallel)
# event = make_event(queue, occupancy)
# inserted = fill(false, length(event))
# scheduled = fill(false, length(event))

# lk = ReentrantLock()
# t = Threads.@spawn while !all(inserted) || !isempty(channel)
#     # println("Inserting")
#     insert_point!(channel, parallel_tree, inserted, lk)
# end
# while !all(inserted)
#     timer = time()
#     points = sum(inserted)
#     println(points)
#     for i in event
#         Threads.@spawn run_event(channel, i, inserted, scheduled, parallel_tree, lk)
#     end
#     # println("Insert per second: ", (sum(inserted)- points)/(time()-timer))
# end