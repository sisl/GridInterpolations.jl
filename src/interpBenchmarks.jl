# This function is used to measure the performance and efficiency of the interpolation functions.

using GridInterpolations


function benchmark(interpType::Type, numDims::Int=6, pointsPerDim::Int=15, numRandomTests::Int=1000; quiet=false, k=1)
    # Calculate interpolation benchmarks
    # create a multi-dimensional test Grid with data, cut points, etc.
	
    # Set up the data structures:
    cutPointMat = 0:(1/(pointsPerDim-1)):1
    cutPoints = [cutPointMat for i=1:numDims]
	
	if interpType == SimplexGrid
		grid = SimplexGrid(tuple(cutPoints...)...)
	elseif interpType == RectangleGrid
		grid = RectangleGrid(tuple(cutPoints...)...)
	elseif interpType == KnnGrid
		grid = KnnGrid(k, tuple(cutPoints...)...)
	end
	
	data = rand(length(grid))

    testPoint = rand(numDims,numRandomTests)

    tic()
    for i=1:numRandomTests
        testM = interpolate(grid,data,testPoint[:,i])
    end
    elapsedTime = toq()

    if !quiet
        println("$(label(grid)) interpolation of $numDims dimensions with $pointsPerDim cut points per dimension:")
        println("  $numRandomTests interpolations required $elapsedTime seconds")
        println(string("  10^6 interpolations would take ", elapsedTime*1000000/numRandomTests, " seconds"))
    end

    return elapsedTime
end


function compareBenchmarks(numDims::Int=6, pointsPerDim::Int=15, numRandomTests::Int=1000; quiet=false,
    returnSpeed=false)
    # Compare rectangular and simplex interpolation benchmarks

    elapsedTimeSimplex = benchmark(SimplexGrid, numDims, pointsPerDim, numRandomTests, quiet=true);
    elapsedTimeRectangle = benchmark(RectangleGrid, numDims, pointsPerDim, numRandomTests, quiet=true);
    elapsedTimeKnn1 = benchmark(KnnGrid, numDims, pointsPerDim, numRandomTests, quiet=true, k=1);
    elapsedTimeKnn3 = benchmark(KnnGrid, numDims, pointsPerDim, numRandomTests, quiet=true, k=3);
#    elapsedTimeKnn5 = benchmark(KnnGrid, numDims, pointsPerDim, numRandomTests, quiet=true, k=5);

    if !quiet
        println("$numRandomTests interpolations of $numDims dimensions with $pointsPerDim cut points per dimension:")
        println("  Rectangle required $elapsedTimeRectangle sec")
        println("  Simplex   required $elapsedTimeSimplex sec")
        println("  1-NN       required $elapsedTimeKnn1 sec")
        println("  3-NN       required $elapsedTimeKnn3 sec")
#        println("  5-NN       required $elapsedTimeKnn5 sec")
    end

    if returnSpeed
        return elapsedTimeRectangle/elapsedTimeSimplex
    end

    return 0.0
end

function compareSpeedUp(nDims, nPoints)
    speed = zeros(nDims)
    for i = 1:nDims
        speed[i] = compareBenchmarks(i, nPoints, 1000, quiet = true, returnSpeed = true)
        gc()
    end
    return speed
end



benchmark(RectangleGrid, quiet=true)
benchmark(RectangleGrid, quiet=true)
benchmark(RectangleGrid)
benchmark(SimplexGrid, quiet=true)
benchmark(SimplexGrid, quiet=true)
benchmark(SimplexGrid)
benchmark(KnnGrid, k=1, quiet=true)
benchmark(KnnGrid, k=1, quiet=true)

Profile.clear()
@profile (for i = 1:100; benchmark(KnnGrid, k=1); end)
Profile.print()
Profile.print(format=:flat)

#compareBenchmarks();
