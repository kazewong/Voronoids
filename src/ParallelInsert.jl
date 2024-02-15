using ThreadsX

function batch_locate(vertices::AbstractArray, tree::DelaunayTree; n_dims::Int = 3)::Vector{Vector{Int}}
    output = Vector{Vector{Int}}(undef, length(vertices))
    for i in 1:length(vertices)
        output[i] = locate(Vector{Int}(), vertices[i], tree; n_dims=n_dims)
    end
    return output
end

function parallel_locate(vertices::Vector{Vector{Float64}}, tree::DelaunayTree; n_dims::Int = 3)::Vector{Vector{Int}}
    batch_size = length(vertices) ÷ Threads.nthreads()
    vertices = Iterators.partition(vertices, batch_size)
    task = map(vertices) do vertice
        Threads.@spawn batch_locate(vertice, tree, n_dims=n_dims)
    end
    return reduce(vcat, fetch.(task))
end

export parallel_locate, batch_locate