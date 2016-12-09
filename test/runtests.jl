using GridInterpolations
using Base.Test
using Grid
using Distances


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
#    showcompact(grid)
#    show(grid)
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



function test_x2ind_knn()
	knn = KnnGrid(3, [2,5,150,10,20,100],[2,5],[1,2,3,5,6])
	@test ind2x(knn, x2ind(knn, [2, 2, 6])) == [2, 2, 6]
	@test ind2x(knn, x2ind(knn, [20, 5, 3])) == [20, 5, 3]
	@test ind2x(knn, x2ind(knn, [2, 2, 1])) == [2, 2, 1]
	@test ind2x(knn, x2ind(knn, [100, 5, 6])) == [100, 5, 6]
end
test_x2ind_knn()


# Check function to get nearest k in dim
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5,10,20,100,150], 12, 2)) == Set([10, 5])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5,10,20,100,150], 11, 2)) == Set([10,5])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5,10,20,100,150], 11, 3)) == Set([10,5,2,20])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5,10,20,100,150], 50, 3)) == Set([5, 10, 20])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5,10,20,100,150], 50, 1)) == Set([20])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5,10,20,100,150], 100, 1)) == Set([100])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5,10,20,100,150], 80, 1)) == Set([100])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5], 1, 1)) == Set([2])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5], 2, 1)) == Set([2])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5], 3, 1)) == Set([2])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5], 4, 1)) == Set([5])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5], 5, 1)) == Set([5])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5], 6, 1)) == Set([5])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5], 3.5, 1)) == Set([2,5])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5], 1, 2)) == Set([2,5])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5], 2, 2)) == Set([2,5])
@test Set(GridInterpolations.get_nearest_k_in_dim([2,5], 3, 2)) == Set([2,5])
@test Set(GridInterpolations.get_nearest_k_in_dim([2], 1, 1)) == Set([2])
@test Set(GridInterpolations.get_nearest_k_in_dim([2], 2, 1)) == Set([2])
@test Set(GridInterpolations.get_nearest_k_in_dim([2], 3, 1)) == Set([2])
@test_throws ErrorException GridInterpolations.get_nearest_k_in_dim([], 3, 1)


# Check implementation of knn
function test_knn_cutpoints_sorted()
    knn = KnnGrid(3, [5,2],[100,5])
    return knn.cutPoints[1] == [2,5] && knn.cutPoints[2] == [5,100]
end
@test test_knn_cutpoints_sorted() == true


# Check that duplicates are disallowed
@test_throws ErrorException KnnGrid(3, [1,1])
@test_throws ErrorException RectangleGrid([1,1])
@test_throws ErrorException SimplexGrid([1,1])


# Sanity check overall implementations
function test_rect_implemented()
    rgrid = RectangleGrid([2,5],[2,5])
	
	@test length(rgrid) == 4
	@test GridInterpolations.dimensions(rgrid) == 2
	
	@test ind2x(rgrid,1) == [2,2]
	@test ind2x(rgrid,2) == [5,2]
	@test ind2x(rgrid,3) == [2,5]
	@test ind2x(rgrid,4) == [5,5]
	x = [0,0]
	ind2x!(rgrid, 4, x)
	@test x == [5,5]
	
	@test interpolants(rgrid, [1,1]) == ([1],[1.0])
	@test interpolants(rgrid, [2,5]) == ([3],[1.0])
	
	indices, weights = interpolants(rgrid, [1.5,3])
	@test indices == [1,3]
	@test isapprox(weights, [0.666667, 0.333333], rtol=1e-5)
	
	@test interpolate(rgrid, [1,2,3,4], [1,1]) == 1.0
	@test maskedInterpolate(rgrid, [1,2,3,4], [1,1], BitArray([false, false, false, false])) == 1.0
	@test isnan(maskedInterpolate(rgrid, [1,2,3,4], [1,1], BitArray([true, false, false, false])))
	
	@test interpolate(rgrid, [1,2,3,4], [1.5,3]) == 1.6666666666666665
	@test maskedInterpolate(rgrid, [1,2,3,4], [1.5,3], BitArray([true, false, false, false])) == 3
	@test maskedInterpolate(rgrid, [1,2,3,4], [1.5,3], BitArray([false, false, true, false])) == 1
end
test_rect_implemented()

function test_knn_implemented()
    knn = KnnGrid(3, [5,2],[100,5])
	@test length(knn) == 4
	@test GridInterpolations.dimensions(knn) == 2
	
	@test ind2x(knn,1) == [2,5]
	@test ind2x(knn,2) == [5,5]
	@test ind2x(knn,3) == [2,100]
	@test ind2x(knn,4) == [5,100]
	x = [0,0]
	ind2x!(knn, 4, x)
	@test x == [5,100]
	
	knn2 = KnnGrid(2, [2,5,150,10,20,100],[2,5],[1,2,3,5,6])
	@test interpolants(knn2, [50,3,3.2]) == ([28,16],[0.5,0.5])
	
    knn3 = KnnGrid(1, [2,5],[2,5])
	@test interpolate(knn3, [1,2,3,4], [2,2]) == 1.0
	@test interpolate(knn3, [1,2,3,4], [5,5]) == 4.0
	@test interpolate(knn3, [1,2,3,4], [6,6]) == 4.0
	@test interpolate(knn3, [1,2,3,4], [3,5]) == 3.0
	
    knn4 = KnnGrid(3, [1,5],[1,5])
	@test isapprox(interpolate(knn4, [1,2,3,4], [2,2]), 2.0, rtol=1e-5)
	@test isapprox(interpolate(knn4, [1,2,3,4], [4,4]), 3.0, rtol=1e-5)
	
	@test maskedInterpolate(knn4, [1,2,3,4], [2,2], BitArray([true, false, false, false])) == 2.5
end
test_knn_implemented()

# Test that we can use a variety of metrics from Distances.jl inside of the fast grid
@test length(KnnFastGrid(3, [2,5], [2,5,10]).balltree.data) == 6
@test length(KnnFastGrid(3, [2,5], [2,5,10], dist_metric=Euclidean()).balltree.data) == 6 # default
@test length(KnnFastGrid(3, [2,5], [2,5,10], dist_metric=Minkowski(3.5)).balltree.data) == 6
@test length(KnnFastGrid(3, [2,5], [2,5,10], dist_metric=Cityblock()).balltree.data) == 6
@test_throws MethodError KnnFastGrid(3, [2,5], [2,5,10], dist_metric=CosineDist()) # cosine distance is a semimetric
@test_throws UndefVarError KnnFastGrid(3, [2,5], [2,5,10], dist_metric=MadeUp())

function test_knnfast_implemented()
    knn = KnnFastGrid(3, [5,2],[100,5])
	@test length(knn) == 4
	@test GridInterpolations.dimensions(knn) == 2
	
	@test ind2x(knn,1) == [5,100]
	@test ind2x(knn,2) == [2,100]
	@test ind2x(knn,3) == [5,5]
	@test ind2x(knn,4) == [2,5]
	x = [0,0]
	ind2x!(knn, 4, x)
	@test x == [2,5]
	
	knn2 = KnnFastGrid(2, [1,2,3,4],[1,2])
	# In this AbstractGrid type, indices are:
	# 1: (1,1), 2: (2,1), 3: (3,1), 4:(4,1)
	# 5: (1,2), 6: (2,2), 7: (3,2), 8:(4,2)
	indices, weights = interpolants(knn2, [1.5, 1])
	@test weights == [0.5,0.5]
	@test Set(indices) == Set([1,2])
	indices, weights = interpolants(knn2, [4, 1.5])
	@test weights == [0.5,0.5]
	@test Set(indices) == Set([4,8])
	
    knn3 = KnnFastGrid(1, [2,5],[2,5])
	# In this AbstractGrid type, indices are:
	# 1: (2,2), 2: (2,5)
	# 3: (5,2), 4: (5,5)
	@test interpolate(knn3, [1,2,3,4], [2,2]) == 1.0
	@test interpolate(knn3, [1,2,3,4], [5,5]) == 4.0
	@test interpolate(knn3, [1,2,3,4], [6,6]) == 4.0
	@test interpolate(knn3, [1,2,3,4], [3,5]) == 3.0
	
    knn4 = KnnFastGrid(3, [1,5],[1,5])
	@test isapprox(interpolate(knn4, [1,2,3,4], [2,2]), 2.0, rtol=1e-5)
	@test isapprox(interpolate(knn4, [1,2,3,4], [4,4]), 3.0, rtol=1e-5)
	
	@test maskedInterpolate(knn4, [1,2,3,4], [2,2], BitArray([true, false, false, false])) == 2.5
end
test_knnfast_implemented()


# check whether rectangle & simplex expect sorted order of cutpoints on dimensions
function test_ordering(grid)	
	@test ind2x(grid,1) == [2,18] # 1
	@test ind2x(grid,2) == [5,18] # 2
	@test ind2x(grid,3) == [2,15] # 3
	@test ind2x(grid,4) == [5,15] # 4
	@test ind2x(grid,5) == [2,12] # 5
	@test ind2x(grid,6) == [5,12] # 6
	
	return interpolate(grid, [1,2,3,4,5,6], [6, 20]) == 2.0
end
@test_throws ErrorException test_ordering( RectangleGrid([2,5], [18,15,12]) )
@test_throws ErrorException test_ordering( SimplexGrid([2,5], [18,15,12]) )
@test_throws ErrorException test_ordering( KnnGrid(k=1[2,5], [18,15,12]) )
@test test_ordering( KnnFastGrid(k=1, [2,5], [18,15,12]) ) == true



#include(joinpath(Pkg.dir("GridInterpolations"), "src", "interpBenchmarks.jl"))
#compareBenchmarks(quiet=true)


println("All tests complete")