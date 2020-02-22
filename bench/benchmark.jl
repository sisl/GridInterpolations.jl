using GridInterpolations
using Statistics
using BenchmarkTools

function benchmark_interpolate(grid, data, rand_points)
    for point in eachcol(rand_points)
        interpolate(grid, data, point)
    end
end

function benchmark(GridType::Type, n_dims::Int=6, n_points_per_dim::Int=15, n_rand_points::Int=100)

    # construct grid
    grid = GridType((range(0, 1, length=n_points_per_dim) for _ in 1:n_dims)...)
    rand_points = rand(n_dims, n_rand_points)
    rand_data = rand(length(grid))

    println("Benchmarking interpolation on $(label(grid))")
    println("dimensionality:          $(n_dims)")
    println("points per dimension:    $(n_points_per_dim)")
    println("number of random points: $(n_rand_points)")

    t = @benchmark benchmark_interpolate($grid, $rand_data, $rand_points) samples=10000
    trial_estimate = median(t)
    println(trial_estimate)

    println("--------------------------------------------------")

    return median(t)
end

for ndims in 1:10

    println("\n##################################################")
    println("### SUMMARY (ndims = $(ndims))")
    println("##################################################")

    trial_estimate_rectangle = benchmark(RectangleGrid, ndims)
    trial_estimate_simplex = benchmark(SimplexGrid, ndims)

    # println(ratio(trial_estimate_rectangle, trial_estimate_simplex))
    println(
        "RectangleGrid vs SimplexGrid: ",
        judge(trial_estimate_rectangle, trial_estimate_simplex)
    )
    println("##################################################\n")

end
