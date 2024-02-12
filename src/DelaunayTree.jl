using AdaptiveKDTrees.KNN
using BoundingSphere

mutable struct DelaunayTree
    vertices::Vector{Vector{Float64}}
    dead::Vector{Bool}
    kdtree::KDTree
    simplices::Vector{Vector{Int}}
    vertices_simplex::Vector{Vector{Int}}
    centers::Vector{Vector{Float64}}
    radii::Vector{Float64}
    neighbors_relation::Vector{Vector{Int}}
end

function circumsphere(vertices::Vector{Vector{Float64}}; n_dims::Int=3)
    if n_dims == 2
        x1, y1 = vertices[1]
        x2, y2 = vertices[2]
        x3, y3 = vertices[3]

        # Midpoints of AB and BC
        D = ((x1 + x2) / 2, (y1 + y2) / 2)
        E = ((x2 + x3) / 2, (y2 + y3) / 2)

        # Slopes of AB and BC
        mAB = (y2 - y1) / (x2 - x1)
        mBC = (y3 - y2) / (x3 - x2)

        # Slopes of perpendicular bisectors
        mD = -1 / mAB
        mE = -1 / mBC

        # Calculating the circumcenter (X, Y)
        X = (mD * D[1] - mE * E[1] + E[2] - D[2]) / (mD - mE)
        Y = mD * (X - D[1]) + D[2]

        # Radius of the circumcircle
        R = sqrt((X - x1)^2 + (Y - y1)^2)

        return ([X, Y], R)
    elseif n_dims == 3
        v1 = vertices[1]
        v2 = vertices[2]
        v3 = vertices[3]
        v4 = vertices[4]

        if (v1==v2) || (v1==v3) || (v1==v4) || (v2==v3) || (v2==v4) || (v3==v4)
            return ((0, 0, 0), 0)
        end

        length_column = [
            v1[1]^2 + v1[2]^2 + v1[3]^2;
            v2[1]^2 + v2[2]^2 + v2[3]^2;
            v3[1]^2 + v3[2]^2 + v3[3]^2;
            v4[1]^2 + v4[2]^2 + v4[3]^2;
        ]

        a = det([v1[1] v1[2] v1[3] 1;
            v2[1] v2[2] v2[3] 1;
            v3[1] v3[2] v3[3] 1;
            v4[1] v4[2] v4[3] 1])

        Dx = det([length_column[1] v1[2] v1[3] 1;
            length_column[2] v2[2] v2[3] 1;
            length_column[3] v3[2] v3[3] 1;
            length_column[4] v4[2] v4[3] 1])

        Dy = - det([length_column[1] v1[1] v1[3] 1;
            length_column[2] v2[1] v2[3] 1;
            length_column[3] v3[1] v3[3] 1;
            length_column[4] v4[1] v4[3] 1])

        Dz = det([length_column[1] v1[1] v1[2] 1;
            length_column[2] v2[1] v2[2] 1;
            length_column[3] v3[1] v3[2] 1;
            length_column[4] v4[1] v4[2] 1])

        center = [Dx/2/a, Dy/2/a, Dz/2/a]
        radius = sqrt((v1[1]-center[1])^2 + (v1[2]-center[2])^2 + (v1[3]-center[3])^2)

        return ([Dx/2/a,Dy/2/a,Dz/2/a], radius) # Return the center coordinates and the radius
    end
end

function initialize_tree_3d(positions::Vector{Vector{Float64}})::DelaunayTree
    center, radius = boundingsphere(positions)
    radius = radius*5
    first_vertex = center + [0, 0, radius]
    second_vertex = center + [radius * cos(0), radius * sin(0), -radius ]
    third_vertex = center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3), -radius]
    fourth_vertex = center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3), -radius]
    ghost_vertex = [center + [0, 0, radius], center + [0, 0, radius], center + [0, 0, radius], center + [radius  * cos(0), radius * sin(0), -radius]]
    vertices = [first_vertex, second_vertex, third_vertex, fourth_vertex, ghost_vertex...]

    simplicies = [[1, 2, 3, 4], [5, 1, 2, 3], [6, 1, 3, 4], [7, 1, 4, 2], [8, 2, 3, 4]]
    vertices_simplex = [[1,2,3,4], [1, 2, 4, 5], [1, 2, 3, 5], [1, 3, 4, 5], [2], [3], [4], [5]]
    centers = [center, [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]
    radii = [radius, 0, 0, 0, 0]
    dead = [false, false, false, false, false]

    neighbors_relation = [[2, 3, 4, 5], [1], [1], [1], [1]]
    kd_tree = KDTree(reduce(hcat, vertices))
    return DelaunayTree(vertices, dead, kd_tree, simplicies, vertices_simplex, centers, radii, neighbors_relation)
end

function initialize_tree_2d(positions::Vector{Vector{Float64}})::DelaunayTree
    center, radius = boundingsphere(positions)
    radius = radius*3
    first_vertex = center + [radius * cos(0* pi / 3), radius * sin(0 * pi / 3)]
    second_vertex = center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3)]
    third_vertex = center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3)]
    ghost_vertex = [center + [radius * cos(0 * pi / 3), radius * sin(0 * pi / 3)], center + [radius * cos(2 * pi / 3), radius * sin(2 * pi / 3)], center + [radius * cos(4 * pi / 3), radius * sin(4 * pi / 3)]]
    vertices = [first_vertex, second_vertex, third_vertex, ghost_vertex...]

    simplicies = [[1, 2, 3], [4, 1, 2], [6, 1, 3], [5, 2, 3]]
    vertices_simplex = [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2], [3], [4]]
    centers = [center, [0, 0], [0, 0], [0, 0]]
    radii = [radius, 0, 0, 0]
    dead = [false, false, false, false]

    neighbors_relation = [[2, 3, 4], [1], [1], [1]]
    kd_tree = KDTree(reduce(hcat, vertices))
    return DelaunayTree(vertices, dead, kd_tree, simplicies, vertices_simplex, centers, radii, neighbors_relation)
end

export initialize_tree_2d, initialize_tree_3d, DelaunayTree, circumsphere