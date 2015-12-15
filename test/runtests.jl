using GridInterpolations
using Base.Test
using Grid


function compareToGrid(testType::Symbol=:random, numDims::Int=3, pointsPerDim::Int=5, testPointsDim::Int=5; numRandomTests::Int=1000, eps::Float64=1e-10)
    # Compare interpolation results from GridInterpolations with the Grid Module:
    # create a multi-dimensional test Grid with data, cut points, etc.
    # (dataM, gridM, gridI) = compareToGrid(2, 3, 5);
    cutsDim = pointsPerDim*ones(Int,numDims)
    gridData = rand(cutsDim...)

    # Set up the data structures:
    gridI = InterpGrid(gridData, BCnearest, InterpLinear)  		# This is for the Grid package
    cutPointMat = 0:(1/(pointsPerDim-1)):1
    cutPoints = [cutPointMat for i=1:numDims]
    gridM = RectangleGrid(tuple(cutPoints...)...) 				# Mykel's interpolation package

    # Package the data in a way that can be read by the interpolation package:
    dataM = Array(Float64,length(gridM))
    for i=1:length(gridM)
        # x = ind2x(gridM, i)
        x = ind2x(gridM, i)
        dataM[i] = gridData[getFractionalIndexes(gridM,x)...]
    end

    # Check a finer grid of coordinates, for comprehensiveness?

    if testType==:random
        # Interpolate in GridInterpolations with the same points, report any discrepancy
        # Select a number of random points and compare the interpolation results.
        testPoint = zeros(numDims)
        for i=1:numRandomTests
            rand!(testPoint)
            testI = gridI[((pointsPerDim-1)*testPoint+1)...]
            testM = interpolate(gridM,dataM,testPoint)

            if (abs(testI-testM)>eps)
                display("Failed Random Point Interpolation Test")
                display("Non-matching interpolation points:")
                display("Test point = $testPoint")
                display("Grid Package interpolation result = $testI")
                display("Interpolation Module result = $testM")
                return false
            end

        end
        return true
    end

    if testType==:extrapNeg
        # Test behavior outside the grid, i.e. extrapolation:
        testPoint = zeros(numDims)
        for i=1:numRandomTests
            testPoint = -rand!(testPoint)
            testI = gridI[((pointsPerDim-1)*testPoint+1)...]
            testM = interpolate(gridM,dataM,testPoint)

            if (abs(testI-testM)>eps)
                display("Failed Random Negative Point Extrapolation Test")
                display("Non-matching interpolation points:")
                display("Test point = $testPoint")
                display("Grid Package interpolation result = $testI")
                display("Interpolation Module result = $testM")
                return false
            end

        end
        return true
    end

    if testType==:extrapPos
        for i=1:numRandomTests
            testPoint = rand(numDims)+1
            testI = gridI[((pointsPerDim-1)*testPoint+1)...]
            testM = interpolate(gridM,dataM,testPoint)

            if (abs(testI-testM)>eps)
                display("Failed Random Positive Point Extrapolation Test")
                display("Non-matching interpolation points:")
                display("Test point = $testPoint")
                display("Grid Package interpolation result = $testI")
                display("Interpolation Module result = $testM")
                return false
            end

        end
        return true
    end

    return true

end

function getFractionalIndexes(g::AbstractGrid, s::Array)
    # Returns the fractional index of sprime within the grid-defined discretization.

    fracInd = Array(Int64,length(g.cutPoints))

    for i=1:length(g.cutPoints)
        gridDisc = g.cutPoints[i]
        fracInd[i] = getFracIndex(gridDisc, s[i])
    end

    return fracInd

end

function getFracIndex(vararray::Array, value::Float64)
    # Searches through the array vararray and returns the fractional index position
    # that value lies in.  It CAPassumes vararray increases monotonically, and returns
    # the first or last element if the value is outside the bounds of the array.

    if value <= vararray[1]
        return 1
    end

    if value >= vararray[end]
        return indmax(vararray)
    end

    i=1
    while value > vararray[i+1]
        i+=1
    end

    return i+(value-vararray[i])÷(vararray[i+1] - vararray[i])

end

function simplexMagic(NDISC::Int=20, NPOINTS::Int=3, checkFileName::AbstractString="", eps::Float64=1e-10)
    if isempty(checkFileName)
        checkFileName = joinpath(Pkg.dir("GridInterpolations"), "test", "simplexMagicTest20.txt")
    end

    val = transpose([ [8.,1,6] [3,5,7] [4,9,2] ]) # transposed magic(3) from matlab

    sGrid = SimplexGrid(1.:NPOINTS, 1.:NPOINTS)
    sInterpVal = zeros(NDISC, NDISC)

    for i = 1:NDISC, j = 1:NDISC
        valX = mapPt(i, NPOINTS, NDISC)
        valY = mapPt(j, NPOINTS, NDISC)
        sInterpVal[i, j] = interpolate(sGrid, vcat(val...), [valX, valY])
    end

    sInterpValTest = readdlm(checkFileName)

    testErr = sum(abs(sInterpVal-sInterpValTest))

    if (testErr > eps)
        display("Failed Simplex Comparison Test")
        display("Test error = $testErr")
    end

    return testErr<eps

end

# constructs a new grid from repr output and tests to see if it's the same
function reprConstruct()
    g1 = SimplexGrid([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
    g2 = eval(parse(repr(g1)))
    return g1.cutPoints == g2.cutPoints
end

function mapPt(pt::Int, NPOINTS::Int, NDISC::Int)
    return pt * (NPOINTS-1) / NDISC + 1
end

function checkCounters()
    dims = 2
    nPts = 100
    grid = RectangleGrid(1.:10, 1.:10)
    if length(grid) != nPts || dimensions(grid) != dims
        return false
    end
    data = rand(10,10)
    temp = interpolate(grid, data, [-1.5,1.5])
    temp = interpolate(grid, data, [1.,1.])
    grid = SimplexGrid(1.:10, 1.:10)
    if length(grid) != nPts || dimensions(grid) != dims
        return false
    end
    temp = interpolate(grid, data, [-1.5,1.5])
    temp = interpolate(grid, data, [1.,1.])
    temp = interpolate(grid, data, [1.5,1.5])
    showcompact(grid)
    show(grid)
    return true
end

function testMask()
    x = [1.5, 1.5]
    vals = collect(1.:9)
    mask = falses(9)
    grid = RectangleGrid(1.:3, 1.:3)
    interped = interpolate(grid, vals, x)
    masked   = maskedInterpolate(grid, vals, x, mask)
    @test_approx_eq_eps interped masked 1e-3
    mask[1] = true # mask the first entry
    truth  = 3.6666
    masked = maskedInterpolate(grid, vals, x, mask)
    @test_approx_eq_eps truth masked 1e-3
    return true
end


@test compareToGrid(:random) == true
@test compareToGrid(:extrapPos) == true
#@test compareToGrid(:extrapNeg)==true

@test simplexMagic() == true

@test reprConstruct() == true

@test checkCounters() == true

@test testMask() == true

include(joinpath(Pkg.dir("GridInterpolations"), "src", "interpBenchmarks.jl"))
rectangleBenchmark(quiet=true)
simplexBenchmark(quiet=true)
compareBenchmarks(quiet=true)
