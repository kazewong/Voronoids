function parallel_locate(vertices::Vector{Vector{Float64}}, tree::DelaunayTree; n_dims::Int = 3)::Vector{Int}
    output = Vector{Vector{Int}}(undef, length(vertices))
    Threads.@threads for i in 1:length(vertices)
        output[i] = locate(Vector{Int}(), vertices[i], tree, n_dims=n_dims)
    end

end