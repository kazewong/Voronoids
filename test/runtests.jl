using Voronoids
using Test
using TimerOutputs
using Random

const tmr = TimerOutput()

function test_2d(n::Int; seed::Int)
    Random.seed!(seed)
    n_dims = 2
    test_points = [rand(n_dims) for i in 1:n]
    tree = initialize_tree_2d(test_points)

    for point in test_points
        insert_point(tree, point, n_dims=n_dims)
    end
    check_delaunay(tree, n_dims=2)

    return tree
end


function test_3d(n::Int; seed::Int)
    Random.seed!(seed)
    n_dims = 3
    test_points = [rand(n_dims) for i in 1:n]
    tree = initialize_tree_3d(test_points)

    for point in test_points
        insert_point(tree, point, n_dims=n_dims)
    end
    check_delaunay(tree, n_dims=3)

    return tree
end

@testset "delaunay test" begin
    @testset "2d" begin
        tree = test_2d(100; seed=123)
    end

    @testset "3d" begin
        tree = test_3d(100; seed=123)
    end
end