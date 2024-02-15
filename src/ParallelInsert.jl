function batch_locate(vertices::AbstractArray, tree::DelaunayTree)::Vector{Vector{Int}}
    output = Vector{Vector{Int}}(undef, length(vertices))
    for i in 1:length(vertices)
        output[i] = locate(Vector{Int}(), vertices[i], tree)
    end
    return output
end

function parallel_locate(vertices::Vector{Vector{Float64}}, tree::DelaunayTree)::Vector{Vector{Int}}
    batch_size = length(vertices) ÷ Threads.nthreads()
    vertices = Iterators.partition(vertices, batch_size)
    task = map(vertices) do vertice
        Threads.@spawn batch_locate(vertice, tree)
    end
    return reduce(vcat, fetch.(task))
end

function identify_nonconflict_points(vertices::Vector{Vector{Float64}}, tree::DelaunayTree)::Vector{Int}
    site_list = parallel_locate(vertices, tree)
    length_site_list = length(site_list)
    output = Vector{Bool}(undef, length_site_list)
    Threads.@threads for i in 1:length_site_list-1
        # output[i] = all(isempty.(map(x->intersect(reduce(vcat,tree.neighbors_relation[site_list[i]]), reduce(vcat, tree.neighbors_relation[site_list[x]])), i+1:length_site_list)))
        output[i] = all(isempty.(map(x->intersect(site_list[i], site_list[x]), i+1:length_site_list)))
    end
    return output
end

export parallel_locate, batch_locate, identify_nonconflict_points