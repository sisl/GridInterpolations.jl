# This function is used to measure the performance and efficiency of the interpolation functions.

using GridInterpolations


function benchmark(interpType::Type, numDims::Int=6, pointsPerDim::Int=15, numRandomTests::Int=1000; quiet=false, k=1)
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
	elseif interpType == KnnGrid
		grid = KnnGrid(k, tuple(cutPoints...)...)
	elseif interpType == KnnFastGrid
		grid = KnnFastGrid(k, tuple(cutPoints...)...)
	end
	
	data = rand(length(grid))

    testPoint = rand(numDims,numRandomTests)

	# Warmup
    interpolate(grid,data,testPoint[:,1])
    tic()
    for i=1:numRandomTests
        testM = interpolate(grid,data,testPoint[:,i])
    end
    elapsedTime = toq()
	@show numDims,pointsPerDim, length(grid), elapsedTime

    if !quiet
        println("$(label(grid)) interpolation of $numDims dimensions with $pointsPerDim cut points per dimension:")
        println("  $numRandomTests interpolations required $elapsedTime seconds")
    end

    return elapsedTime
end


function compareBenchmarks(numDims::Int=6, pointsPerDim::Int=15, numRandomTests::Int=1000; quiet=false)
    # Compare rectangular and simplex interpolation benchmarks

    elapsedTimeSimplex = benchmark(SimplexGrid, numDims, pointsPerDim, numRandomTests, quiet=true);
    elapsedTimeRectangle = benchmark(RectangleGrid, numDims, pointsPerDim, numRandomTests, quiet=true);
    elapsedTimeKnn1 = benchmark(KnnGrid, numDims, pointsPerDim, numRandomTests, quiet=true, k=1);
    elapsedTimeKnn3 = benchmark(KnnGrid, numDims, pointsPerDim, numRandomTests, quiet=true, k=3);
#    elapsedTimeKnn5 = benchmark(KnnGrid, numDims, pointsPerDim, numRandomTests, quiet=true, k=5);
    elapsedTimeKnnFast1 = benchmark(KnnFastGrid, numDims, pointsPerDim, numRandomTests, quiet=true, k=1);
    elapsedTimeKnnFast3 = benchmark(KnnFastGrid, numDims, pointsPerDim, numRandomTests, quiet=true, k=3);

    if !quiet
        println("$numRandomTests interpolations of $numDims dimensions with $pointsPerDim cut points per dimension:")
        println("  Rectangle required $elapsedTimeRectangle sec")
        println("  Simplex   required $elapsedTimeSimplex sec")
        println("  1-NN(m)   required $elapsedTimeKnn1 sec")
        println("  3-NN(m)   required $elapsedTimeKnn3 sec")
#        println("  5-NN       required $elapsedTimeKnn5 sec")
        println("  1-NN(b)   required $elapsedTimeKnnFast1 sec")
        println("  3-NN(b)   required $elapsedTimeKnnFast3 sec")
    end

end

function compareSpeedUp(nDims, nPoints; marginoferror=0.1)
	# Holds the number of nPoints constant to within 
	# marginoferror, varying the number of dimensions,
	# and returns the speedup of simplex over the
	# rectangle interpolation
	
	nTests = 1000
	
	# Compile everything & warm up the cache
	benchmark(SimplexGrid, 1, 10, nTests, quiet=true);
	benchmark(SimplexGrid, 1, 10, nTests, quiet=true);
	benchmark(RectangleGrid, 1, 10, nTests, quiet=true);
	benchmark(RectangleGrid, 1, 10, nTests, quiet=true);

	println("begin")
    speedup = []
    for i = 1:nDims
		nPointsPerDim = nPoints^(1/i)
		nPointsPerDim = convert(Int, round(nPointsPerDim))
		
		if abs(nPointsPerDim^i - nPoints) > marginoferror*nPoints
			continue
		end
		
		sspeed = benchmark(SimplexGrid, i, nPointsPerDim, nTests, quiet=true);
		rspeed = benchmark(RectangleGrid, i, nPointsPerDim, nTests, quiet=true);
		
		@show sspeed, rspeed, rspeed/sspeed
		push!(speedup, (i,rspeed/sspeed) )
        gc()
    end
    return speedup
end


compareBenchmarks(8, 10);
#compareSpeedUp(1000,100000000)

#=

# Warm the cache & get results
benchmark(RectangleGrid, quiet=true)
benchmark(RectangleGrid, quiet=true)
benchmark(RectangleGrid)

# Warm the cache & get results
benchmark(SimplexGrid, quiet=true)
benchmark(SimplexGrid, quiet=true)
benchmark(SimplexGrid)

# Warm the cache & get results
benchmark(KnnGrid, k=1, quiet=true)
benchmark(KnnGrid, k=1, quiet=true)
Profile.clear()
@profile (for i = 1:100; benchmark(KnnGrid, k=1); end)
Profile.print()
Profile.print(format=:flat)

=#