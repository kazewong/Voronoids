using LinearAlgebra
struct Point
    pos::Vector{Float64}
end

function solve_plane(points::Vector{Point})::Vector{Float64}
    distance = transpose(reduce(hcat,map((x) -> x.pos .- points[1].pos, points[2:end])))
    coefs = Array{Float64}(undef, length(points) + 1)
    for i in 1:length(points)
        unit_vector = zeros(length(points))
        unit_vector[i] = 1
        coefs[i] = det(vcat(distance, reshape(unit_vector, (1, size(unit_vector, 1)))))
    end
    coefs[end] = -dot(coefs[1:end-1], points[1].pos)
    return coefs
end

function distance_from_plane(test_point::Point, plane::Vector{Point})::Float64
    plane_coefs = solve_plane(plane)
    return abs(dot(plane_coefs[1:end-1], test_point.pos) + plane_coefs[end]) / norm(plane_coefs[1:end-1])
end

points = Vector{Point}(undef, 10)
for i in 1:10
    points[i] = Point(rand(4))
end

solve_plane(points[1:4])
dist = distance_from_plane(points[5], points[1:4])

function FindHull(points::Array{Point,1})
    # ...
end