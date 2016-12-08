# This function is used to measure the performance and efficiency of the interpolation functions.

using GridInterpolations


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
	elseif interpType == KnnGrid
		grid = KnnGrid(k, tuple(cutPoints...)...)
	elseif interpType == KnnFastGrid
		grid = KnnFastGrid(k, tuple(cutPoints...)...)
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
		tic()
		for i=1:totalInterpolations
			interpolate(grid,data,testPoint[:,i])
		end
		elapsedTime = toq()
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
    elapsedTimeKnn1, sdKnn1 = benchmark(KnnGrid, numDims, pointsPerDim, numRandomTests, quiet=true, k=1);
    elapsedTimeKnn3, sdKnn3 = benchmark(KnnGrid, numDims, pointsPerDim, numRandomTests, quiet=true, k=3);
#    elapsedTimeKnn5 = benchmark(KnnGrid, numDims, pointsPerDim, numRandomTests, quiet=true, k=5);
    elapsedTimeKnnFast1, sdKnnFast1 = benchmark(KnnFastGrid, numDims, pointsPerDim, numRandomTests, quiet=true, k=1);
    elapsedTimeKnnFast3, sdKnnFast3 = benchmark(KnnFastGrid, numDims, pointsPerDim, numRandomTests, quiet=true, k=3);

    if !quiet
        println("$numRandomTests interpolations of $numDims dimensions with $pointsPerDim cut points per dimension:")
        println("  Rectangle required $elapsedTimeRectangle +/- $sdRectangle sec")
        println("  Simplex   required $elapsedTimeSimplex +/- $sdSimplex sec")
        println("  1-NN(m)   required $elapsedTimeKnn1 +/- $sdKnn1 sec")
        println("  3-NN(m)   required $elapsedTimeKnn3 +/- $sdKnn3 sec")
#        println("  5-NN       required $elapsedTimeKnn5 sec")
        println("  1-NN(b)   required $elapsedTimeKnnFast1 +/- $sdKnnFast1 sec")
        println("  3-NN(b)   required $elapsedTimeKnnFast3 +/- $sdKnnFast3 sec")
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
		
		sspeed, ssd = benchmark(SimplexGrid, i, nPointsPerDim, nTests, quiet=true);
		rspeed, rsd = benchmark(RectangleGrid, i, nPointsPerDim, nTests, quiet=true);
		
		@show sspeed, rspeed, rspeed/sspeed
		push!(speedup, (i,rspeed/sspeed) )
        gc()
    end
    return speedup
end


#compareBenchmarks(8, 10);
#compareSpeedUp(1000,100000000)

numPts = 10
for numDim = 1:10
	try
		mean, std = benchmark(RectangleGrid, numDim, numPts, 5000, quiet=true)
		println("rectangle $numDim $mean $std ")
	catch
		# do nothing
	end
	
	try
		mean, std = benchmark(SimplexGrid, numDim, numPts, 5000, quiet=true)
		println("simplex $numDim $mean $std ")
	catch
		# do nothing
	end
	
	try
		mean, std = benchmark(KnnGrid, numDim, numPts, 5000, quiet=true, k=1)
		println("knn $numDim $mean $std ")
	catch
		# do nothing
	end
	
	try
		mean, std = benchmark(KnnGrid, numDim, numPts, 5000, quiet=true, k=1)
		println("knn $numDim $mean $std ")
	catch
		# do nothing
	end
	
	try
		mean, std = benchmark(KnnFastGrid, numDim, numPts, 5000, quiet=true, k=3)
		println("knnfast $numDim $mean $std ")
	catch
		# do nothing
	end
	
	try
		mean, std = benchmark(KnnFastGrid, numDim, numPts, 5000, quiet=true, k=3)
		println("knnfast $numDim $mean $std ")
	catch
		# do nothing
	end
	
end



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
benchmark(KnnGrid, quiet=true)
benchmark(KnnGrid, quiet=true)
benchmark(KnnGrid)

# Warm the cache & get results
benchmark(KnnFastGrid, quiet=true)
benchmark(KnnFastGrid, quiet=true)
benchmark(KnnFastGrid)

# Warm the cache & get results
benchmark(KnnGrid, k=1, quiet=true)
benchmark(KnnGrid, k=1, quiet=true)
Profile.clear()
@profile (for i = 1:100; benchmark(KnnGrid, k=1); end)
Profile.print()
Profile.print(format=:flat)

=#