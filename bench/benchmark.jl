using GridInterpolations
using Statistics
using BenchmarkTools

function benchmark_interpolate(grid, data, rand_points)
    for point in eachcol(rand_points)
        interpolate(grid, data, point)
    end
end

function benchmark(GridType::Type, n_dims::Int; n_points_per_dim::Int=15, n_rand_points::Int=100, verbosity::Int=1)

    # construct grid
    grid = GridType((range(0, 1, length=n_points_per_dim) for _ in 1:n_dims)...)
    rand_points = rand(n_dims, n_rand_points)
    rand_data = rand(length(grid))

    trial = @benchmark benchmark_interpolate($grid, $rand_data, $rand_points) samples=10000

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

verbosity = 1
for ndims in 1:6

    println("\n##################################################")
    println("### Benchmark for ndims = $(ndims)")
    println("##################################################")

    trial_rectangle = benchmark(RectangleGrid, ndims; verbosity=verbosity)
    trial_simplex = benchmark(SimplexGrid, ndims; verbosity=verbosity)

    # println(ratio(trial_estimate_rectangle, trial_estimate_simplex))
    println(
        "RectangleGrid vs SimplexGrid: ",
        judge(median(trial_rectangle), median(trial_simplex))
    )
    println("##################################################\n")

end
