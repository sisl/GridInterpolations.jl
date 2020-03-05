using GridInterpolations
using Statistics
using BenchmarkTools

function _benchmark_interpolate(grid, data, rand_points)
    for point in eachcol(rand_points)
        interpolate(grid, data, point)
    end
end

function benchmark_interpolate(GridType::Type, n_dims::Int; n_points_per_dim::Int=15, n_rand_points::Int=100, verbosity::Int=1)

    # construct grid
    grid = GridType((range(0, 1, length=n_points_per_dim) for _ in 1:n_dims)...)
    rand_points = rand(n_dims, n_rand_points)
    rand_data = rand(length(grid))

    trial = @benchmark _benchmark_interpolate($grid, $rand_data, $rand_points) samples=10000

    if verbosity > 0
        println("Benchmarking interpolation on $(label(grid))")
        println("dimensionality:          $(n_dims)")
        println("points per dimension:    $(n_points_per_dim)")
        println("number of random points: $(n_rand_points)")
        if verbosity > 1
            show(stdout, MIME"text/plain"(), trial)  # make sure printing is verbose by using MIME
        else
            println("\n", trial)
        end
        println("\n--------------------------------------------------")
    end

    return trial
end

function _benchmark_ind2x(grid, indices)
    x = zeros(ndims(grid))
    for i in indices
        x += ind2x(grid, i)
    end
end

function benchmark_ind2x(n_dims=6, n_points_per_dim=10)
    grid = RectangleGrid((range(0, 1, length=n_points_per_dim) for _ in 1:n_dims)...)
    trial = @benchmark _benchmark_ind2x($grid, $(1:length(grid))) samples=1000
end

function benchmark(verbosity = 1)

    for ndims in 1:6

        println("\n##################################################")
        println("### Benchmark for ndims = $(ndims)")
        println("##################################################")

        trial_rectangle = benchmark_interpolate(RectangleGrid, ndims; verbosity=verbosity)
        trial_simplex = benchmark_interpolate(SimplexGrid, ndims; verbosity=verbosity)

        # println(ratio(trial_estimate_rectangle, trial_estimate_simplex))
        println(
            "RectangleGrid vs SimplexGrid: ",
            judge(median(trial_rectangle), median(trial_simplex))
        )

        println("ind2x:", benchmark_ind2x())
        println("##################################################\n")


    end

end

benchmark()
