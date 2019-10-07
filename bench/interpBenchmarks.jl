# This function is used to measure the performance and efficiency of the interpolation functions.

using GridInterpolations
using Statistics

function benchmark(interpType::Type, numDims::Int=6, pointsPerDim::Int=15, numRandomTests::Int=5000; quiet=false, k=1)
    # Calculate interpolation benchmarks
    # create a multi-dimensional test Grid with data, cut points, etc.
    
    # Set up the data structures:
    cutPointMat = 0:(1/(pointsPerDim-1)):1
    cutPoints = [cutPointMat for i=1:numDims]
    cut_counts = Int[length(cutPoints[i]) for i = 1:length(cutPoints)]
  
    if interpType == SimplexGrid
        grid = SimplexGrid(tuple(cutPoints...)...)
    elseif interpType == RectangleGrid
        grid = RectangleGrid(tuple(cutPoints...)...)
    end
  
    data = rand(length(grid))


    # Warmup
    interpolate(grid,data,rand(numDims))
  
    # Run the evaluation: totalInterpolations attempts
    # at interpolating numRandomTests points
    totalInterpolations = 1000
  
    timesToInterpolate = zeros(numRandomTests)
    
    averageTime = 0
    for test=1:numRandomTests
    testPoint = rand(numDims, totalInterpolations)
    elapsedTime = @elapsed begin
        for i=1:totalInterpolations
            interpolate(grid,data,testPoint[:,i])
        end
    end
    timesToInterpolate[test] = elapsedTime
    end
    averageTime = mean(timesToInterpolate)
    plusMinus95CI = std(timesToInterpolate)*2
    #@show numDims,pointsPerDim, length(grid), averageTime

    if !quiet
        println("$(label(grid)) interpolation of $numDims dimensions with $pointsPerDim cut points per dimension:")
        println("  requires $averageTime seconds (+/- $plusMinus95CI for 95% CI) on average to interpolate $totalInterpolations test points (averaged over $numRandomTests attempts)")
    end

    return averageTime, plusMinus95CI
end


function compareBenchmarks(numDims::Int=6, pointsPerDim::Int=15, numRandomTests::Int=1000; quiet=false)
    # Compare rectangular and simplex interpolation benchmarks

    elapsedTimeSimplex, sdSimplex = benchmark(SimplexGrid, numDims, pointsPerDim, numRandomTests, quiet=true);
    elapsedTimeRectangle, sdRectangle = benchmark(RectangleGrid, numDims, pointsPerDim, numRandomTests, quiet=true);

    println("$numRandomTests interpolations of $numDims dimensions with $pointsPerDim cut points per dimension:")
    println("  Rectangle required $elapsedTimeRectangle +/- $sdRectangle sec")
    println("  Simplex   required $elapsedTimeSimplex +/- $sdSimplex sec")

end

function compareSpeedUp(nDims, nPoints; marginoferror=0.1)
    # Holds the number of nPoints constant to within 
    # marginoferror, varying the number of dimensions,
    # and returns the speedup of simplex over the
    # rectangle interpolation
  
    nTests = 30
  
    # Compile everything
    benchmark(RectangleGrid, 3, 10, nTests, quiet=true);

    speedup = []
    println("How large is the simplex grid speed up over the multilinear grid?")
    for i = 2:nDims
        nPointsPerDim = nPoints^(1/i)
        nPointsPerDim = convert(Int, round(nPointsPerDim))
    
        if abs(nPointsPerDim^i - nPoints) > marginoferror*nPoints
            continue
        end
    
        sspeed, ssd = benchmark(SimplexGrid, i, nPointsPerDim, nTests, quiet=true);
        println("  limiting to $i dimensions and therefore $nPointsPerDim points per dim:")
        println("    mean speed: $sspeed, std dev: $ssd")

    
        #push!(speedup, (i,rspeed/sspeed) )
        GC.gc()
    end
    return speedup
end





compareBenchmarks();
compareSpeedUp(1000,100000)

# Warm the cache & get results
#benchmark(RectangleGrid, quiet=true)
#benchmark(RectangleGrid, quiet=true)
#benchmark(RectangleGrid)

# Warm the cache & get results
#benchmark(SimplexGrid, quiet=true)
#benchmark(SimplexGrid, quiet=true)
#benchmark(SimplexGrid)

#=
Profile.clear()
@profile (for i = 1:100; benchmark(SimplexGrid, k=1); end)
Profile.print()
Profile.print(format=:flat)
=#
