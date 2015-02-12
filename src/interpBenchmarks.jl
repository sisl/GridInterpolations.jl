# This function is used to measure the performance and efficiency of the interpolation functions.

using GridInterpolations



function simplexBenchmark(numDims::Int=6, pointsPerDim::Int=15, numRandomTests::Int=1000; quiet=false)
    # Calculate simplex interpolation benchmarks
    # create a multi-dimensional test Grid with data, cut points, etc.

    # cutsDim = pointsPerDim*ones(Int,numDims)
    # gridData = rand(cutsDim...)

    # Set up the data structures:
    cutPointMat = [0:(1/(pointsPerDim-1)):1]
    cutPoints = [cutPointMat for i=1:numDims]
    sGrid = SimplexGrid(tuple(cutPoints...)...)
    sData = rand(length(sGrid))

    testPoint = rand(numDims,numRandomTests)

    tic()
    for i=1:numRandomTests
        testM = interpolate(sGrid,sData,testPoint[:,i])
    end
    elapsedTime = toc()

    if !quiet
        display("Simplex interpolation of $numDims dimensions with $pointsPerDim cut points per dimenion:")
        display("$numRandomTests interpolations required $elapsedTime seconds")
        display(string("10^6 interpolations would take ", elapsedTime*1000000/numRandomTests, " seconds"))
    end

    return elapsedTime

end


function rectangleBenchmark(numDims::Int=6, pointsPerDim::Int=15, numRandomTests::Int=1000; quiet=false)
    # Calculate rectangular interpolation benchmarks
    # create a multi-dimensional test Grid with data, cut points, etc.

    # May want to include the Grid module as well, to compare with this implementation

    # Set up the data structures:
    cutPointMat = [0:(1/(pointsPerDim-1)):1]
    cutPoints = [cutPointMat for i=1:numDims]
    rGrid = RectangleGrid(tuple(cutPoints...)...)
    rData = rand(length(rGrid))

    testPoint = rand(numDims,numRandomTests)

    tic()
    for i=1:numRandomTests
        testM = interpolate(rGrid,rData,testPoint[:,i])
    end
    elapsedTime = toc()

    if !quiet
        display("Rectangular interpolation of $numDims dimensions with $pointsPerDim cut points per dimenion:")
        display("$numRandomTests interpolations required $elapsedTime seconds")
        display(string("10^6 interpolations would take ", elapsedTime*1000000/numRandomTests, " seconds"))
    end
    return elapsedTime

end

function compareBenchmarks(numDims::Int=6, pointsPerDim::Int=15, numRandomTests::Int=1000; quiet=false)
    # Compare rectangular and simplex interpolation benchmarks

    cutsDim = pointsPerDim*ones(Int,numDims)
    gridData = rand(cutsDim...)

    # Set up the data structures:
    cutPointMat = [0:(1/(pointsPerDim-1)):1]
    cutPoints = [cutPointMat for i=1:numDims]
    rGrid = RectangleGrid(tuple(cutPoints...)...)
    sGrid = SimplexGrid(tuple(cutPoints...)...)
    rsData = rand(length(sGrid))

    testPoint = rand(numDims,numRandomTests)

    tic()
    for i=1:numRandomTests
        testM = interpolate(sGrid,rsData,testPoint[:,i])
    end
    elapsedTimeSimplex = toc()

    tic()
    for i=1:numRandomTests
        testM = interpolate(rGrid,rsData,testPoint[:,i])
    end
    elapsedTimeRectangle = toc()

    if !quiet
        display("$numRandomTests interpolations of $numDims dimensions with $pointsPerDim cut points per dimenion:")
        display("Simplex   required $elapsedTimeSimplex sec")
        display("Rectangle required $elapsedTimeRectangle sec")
        display(string("Simplex was faster by a factor of ", elapsedTimeRectangle/elapsedTimeSimplex))
    end

    return 0

end
