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

test_points = [rand(n_dims) for i in 1:n]
test_points2 = [rand(n_dims) for i in 1:1e5]
if n_dims == 3
    tree = initialize_tree_3d(test_points)
else
    tree = initialize_tree_2d(test_points)
end
sites = locate(Vector{Int}(), test_points2[1], tree)

for i in 1:n
    add_vertex!(tree, test_points[i], n_dims=n_dims)
end

parallel_tree = deepcopy(tree)

site_list, neighbor_list, occupancy = identify_conflicts(test_points2[1:128], tree)
groups = group_points(site_list, neighbor_list, occupancy)

const update_channel = Channel{TreeUpdate}(128);

function schedule()
    for group in groups
        if length(group) == 1
            update = make_update(test_points[group[1]], site_list[group[1]], tree, n_dims=n_dims)
            put!(update_channel, update)
        end
    end
end

@async schedule()

b = map(x->neighbor_list[x[1]],filter(x->length(x)==1, groups))
c = map(x->site_list[x[1]],filter(x->length(x)==1, groups))
all(map(y->all(isempty.(map(x->intersect(b[y],b[x]), deleteat!(collect(1:length(b)),y)))), 1:length(b)))
all(map(y->all(isempty.(map(x->intersect(c[y],c[x]), deleteat!(collect(1:length(c)),y)))), 1:length(c)))
all(map(y->all(isempty.(map(x->intersect(b[y],c[x]), deleteat!(collect(1:length(c)),y)))), 1:length(b)))

killed = Vector{Int}()
for i in 1:update_channel.n_avail_items
    update = take!(update_channel)
    # killed = vcat(killed, update.killed_sites)
    insert_point!(tree, update)
end