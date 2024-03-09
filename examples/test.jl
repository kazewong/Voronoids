using Revise
using Voronoids
using TimerOutputs
using BenchmarkTools
using ThreadPinning
using Random

pinthreads(:cores)

const tmr = TimerOutput()

Random.seed!(1234)

n = 1000
n_dims = 3
batch_size = 1000

test_points = [rand(n_dims) for i in 1:n]
test_points2 = [rand(n_dims) for i in 1:1e4]
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

test_points3 = [rand(n_dims) for i in 1:1e4]
inserted = fill(false, length(test_points3))

n_parallel = 10000

occupancy = Dict{Int, Vector{Int}}()
queue = make_queue!(test_points3, occupancy, parallel_tree)
channel = Channel{Tuple{Int,TreeUpdate}}(n_parallel)
event = make_event(queue, occupancy)
inserted = fill(false, length(event))
scheduled = fill(false, length(event))

lk = ReentrantLock()
t = @async while !all(inserted) || !isempty(channel)
    # println("Inserting")
    insert_point!(channel, parallel_tree, inserted, lk)
end
while !all(inserted)
    timer = time()
    points = sum(inserted)
    queue_index = event[inserted.==false]
    queue_index = event[1:min(length(queue_index), n_parallel)]
    for i in queue_index
        run_event(channel, i, inserted, scheduled, parallel_tree)
    end
    println("Insert per second: ", (sum(inserted)- points)/(time()-timer))
end