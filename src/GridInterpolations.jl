module GridInterpolations

using Distances # for distance metrics for knn
using Iterators # for product function used in knn

export AbstractGrid, RectangleGrid, SimplexGrid, KnnGrid, dimensions, length, label, ind2x, ind2x!, x2ind, interpolate, maskedInterpolate, interpolants

abstract AbstractGrid

type RectangleGrid <: AbstractGrid
    cutPoints::Vector{Vector{Float64}}
    cut_counts::Vector{Int}
    cuts::Vector{Float64}
    index::Vector{Int}
    weight::Vector{Float64}
    index2::Vector{Int}
    weight2::Vector{Float64}

    function RectangleGrid(cutPoints...)
        cut_counts = Int[length(cutPoints[i]) for i = 1:length(cutPoints)]
        cuts = vcat(cutPoints...)
        myCutPoints = Array(Vector{Float64}, length(cutPoints))
        for i = 1:length(cutPoints)
			if length(Set(cutPoints[i])) != length(cutPoints[i])
				error(@sprintf("Duplicates cutpoints are not allowed (duplicates observed in dimension %d)",i))
			end
            myCutPoints[i] = cutPoints[i]
        end
        numDims = length(cutPoints)
        index = zeros(Int, 2^numDims)
        weight = zeros(Float64, 2^numDims)
        index[1] = 1
        weight[1] = 1.0
        index2 = zeros(Int, 2^numDims)
        weight2 = zeros(Float64, 2^numDims)
        index2[1] = 1
        weight2[1] = 1.0
        new(myCutPoints, cut_counts, cuts, index, weight, index2, weight2)
    end
end

type SimplexGrid <: AbstractGrid
    cutPoints::Vector{Vector{Float64}}
    cut_counts::Vector{Int}
    cuts::Vector{Float64}
    index::Vector{Int}
    weight::Vector{Float64}
    x_p::Vector{Float64} # residuals
    ihi::Vector{Int} # indices of cuts above point
    ilo::Vector{Int} # indices of cuts below point
    n_ind::Vector{Int}

    function SimplexGrid(cutPoints...)
        cut_counts = Int[length(cutPoints[i]) for i = 1:length(cutPoints)]
        cuts = vcat(cutPoints...)
        myCutPoints = Array(Vector{Float64}, length(cutPoints))
        for i = 1:length(cutPoints)
			if length(Set(cutPoints[i])) != length(cutPoints[i])
				error(@sprintf("Duplicates cutpoints are not allowed (duplicates observed in dimension %d)",i))
			end
            myCutPoints[i] = cutPoints[i]
        end
        numDims = length(cutPoints)
        index = zeros(Int, numDims+1) # d+1 points for simplex
        weight = zeros(Float64, numDims+1)
        x_p = zeros(numDims) # residuals
        ihi = zeros(Int, numDims) # indicies of cuts above point
        ilo = zeros(Int, numDims) # indicies of cuts below point
        n_ind = zeros(Int, numDims)
        new(myCutPoints, cut_counts, cuts, index, weight, x_p, ihi, ilo, n_ind)
    end
end

type KnnGrid <: AbstractGrid
    cutPoints::Vector{Vector{Float64}} # points in each dimension, guaranteed sorted
    cut_counts::Vector{Int} # count of cuts in each dimension
#    index::Vector{Int}
#    weight::Vector{Float64}
    k::Int

    function KnnGrid(k, cutPoints...)
        myCutPoints = Array(Vector{Float64}, length(cutPoints))
        cut_counts = Int[length(cutPoints[i]) for i = 1:length(cutPoints)]
        for i = 1:length(cutPoints)
            if length(Set(cutPoints[i])) != length(cutPoints[i])
                error(@sprintf("Duplicates cutpoints are not allowed (duplicates observed in dimension %d)",i))
            end
            myCutPoints[i] = sort(cutPoints[i])
        end
#        numDims = length(cutPoints)
#        index = zeros(Int, 2^numDims)
#        weight = zeros(Float64, 2^numDims)
        
        new(myCutPoints, cut_counts, k)#, index, weight, k)
    end
end
Base.length(grid::RectangleGrid) = prod(grid.cut_counts)
Base.length(grid::SimplexGrid) = prod(grid.cut_counts)
Base.length(grid::KnnGrid) = prod(grid.cut_counts)

dimensions(grid::RectangleGrid) = length(grid.cut_counts)
dimensions(grid::SimplexGrid) = length(grid.cut_counts)
dimensions(grid::KnnGrid) = length(grid.cut_counts)

label(grid::RectangleGrid) = "multilinear interpolation grid"
label(grid::SimplexGrid) = "simplex interpolation grid"
label(grid::KnnGrid) = "knn interpolation grid"

Base.showcompact(io::IO, grid::AbstractGrid) = print(io, "$(typeof(grid)) with $(length(grid)) points")
Base.show(io::IO, grid::AbstractGrid) = Base.showcompact(io, grid)

# showall returns a valid constructor incantation - it will be called when repr is called on a grid
function Base.showall(io::IO, grid::AbstractGrid)
    print(io, "$(typeof(grid))(")
    for v in grid.cutPoints
        show(io, v)
        print(io, ',')
    end
    print(io, ')')
end

function ind2x(grid::AbstractGrid, ind::Int)
    ndims = dimensions(grid)
    x = Array(Float64, ndims)
    ind2x!(grid, ind, x)
    x::Array{Float64}
end

function ind2x!(grid::AbstractGrid, ind::Int, x::Array)
    # Populates x with the value at ind. 
	# In-place version of ind2x.
	# Example:
	#   rgrid = RectangleGrid([2,5],[20,50])
	#   x = [0,0]
	#   ind2x!(rgrid,4,x)  # x now contains [5,50]
    #   @show x            # displays [5,50]
    ndims = dimensions(grid)
    stride = grid.cut_counts[1]
    for i=2:ndims-1
        stride *= grid.cut_counts[i]
    end

    for i=(ndims-1):-1:1
        rest = rem(ind-1, stride) + 1
        x[i + 1] = grid.cutPoints[i + 1][div(ind - rest, stride) + 1]
        ind = rest
        stride = div(stride, grid.cut_counts[i])
    end
    x[1] = grid.cutPoints[1][ind]
    nothing
end

function x2ind(grid::KnnGrid, x::Vector)
    # Inverse of ind2x(grid, idx) = vector
    # Here, we provide a vector x and the grid and get the corresponding index
	# Uses a rectangular grid & relies on the cutPoints being sorted,
	# so only valid for KnnGrid
    num_dim = dimensions(grid)
    totalDimensionsPassed = 1
    index = 1
    for dimid in 1:num_dim
        idx = searchsorted(grid.cutPoints[dimid], x[dimid])
        index += totalDimensionsPassed*(collect(idx) - 1)
        totalDimensionsPassed *= length(grid.cutPoints[dimid])
    end
    index[1]::Int
end


function get_nearest_k_in_dim(this_dimension::Vector, x::Number, k::Int)
    # Assumes this_dimension is sorted
    # Returns the k closest values to x in this_dimension
    # If the k'th closest value is tied for distance to x, will return k+1 values
	if length(this_dimension) == 0
		error(@sprintf("Dimension must have one or more cutpoints"))
	end
	
    possibilities_this_dim = []
    
    firstIndexGte = searchsortedfirst(this_dimension, x)
	
    left = firstIndexGte - 1
    right = firstIndexGte	
    while length(possibilities_this_dim) < k && (left >= 1 || right <= length(this_dimension))
        if left < 1
            push!(possibilities_this_dim, this_dimension[right])
            right += 1
        elseif right > length(this_dimension)
            push!(possibilities_this_dim, this_dimension[left])
            left -= 1
                
        elseif (x - this_dimension[left]) == (this_dimension[right] - x) 
            push!(possibilities_this_dim, this_dimension[left])
            left -= 1
            push!(possibilities_this_dim, this_dimension[right])
            right += 1
        elseif (x - this_dimension[left]) < (this_dimension[right] - x) 
            push!(possibilities_this_dim, this_dimension[left])
            left -= 1
        else
            push!(possibilities_this_dim, this_dimension[right])
            right += 1
        end
    end

    return possibilities_this_dim
end


# masked interpolation ignores points that are masked
function maskedInterpolate(grid::AbstractGrid, data::DenseArray, x::Vector, mask::BitArray{1})
    index, weight = interpolants(grid, x)
    val = 0
    totalWeight = 0
    for i = 1:length(index)
        if mask[index[i]]
            continue
        end
        val += data[index[i]] * weight[i]
        totalWeight += weight[i]
    end
    return val / totalWeight
end

interpolate(grid::AbstractGrid, data::Matrix, x::Vector) = interpolate(grid, map(Float64, data[:]), x)

function interpolate(grid::AbstractGrid, data::DenseArray, x::Vector)
    index, weight = interpolants(grid, x)
    dot(data[index], weight)
end

function interpolants(grid::RectangleGrid, x::Vector)
    cut_counts = grid.cut_counts
    cuts = grid.cuts

    # Reset the values in index and weight:
    fill!(grid.index,0)
    fill!(grid.index2,0)
    fill!(grid.weight,0)
    fill!(grid.weight2,0)
    grid.index[1] = 1
    grid.index2[1] = 1
    grid.weight[1] = 1.
    grid.weight2[1] = 1.

    l = 1
    subblock_size = 1
    cut_i = 1
    n = 1
    for d = 1:length(x)
        coord = x[d]
        lasti = cut_counts[d]+cut_i-1
        ii = cut_i

        if coord <= cuts[ii]
            i_lo, i_hi = ii, ii
        elseif coord >= cuts[lasti]
            i_lo, i_hi = lasti, lasti
        else
            while cuts[ii] < coord
                ii = ii + 1
            end
            if cuts[ii] == coord
                i_lo, i_hi = ii, ii
            else
                i_lo, i_hi = (ii-1), ii
            end
        end

        if i_lo == i_hi
            for i = 1:l
                grid.index[i] += (i_lo - cut_i)*subblock_size
            end
        else
            low = (1 - (coord - cuts[i_lo])/(cuts[i_hi]-cuts[i_lo]))
            for i = 1:l
                grid.index2[i  ] = grid.index[i] + (i_lo-cut_i)*subblock_size
                grid.index2[i+l] = grid.index[i] + (i_hi-cut_i)*subblock_size
            end
            copy!(grid.index,grid.index2)
            for i = 1:l
                grid.weight2[i  ] = grid.weight[i]*low
                grid.weight2[i+l] = grid.weight[i]*(1-low)
            end
            copy!(grid.weight,grid.weight2)
            l = l*2
            n = n*2
        end
        cut_i = cut_i + cut_counts[d]
        subblock_size = subblock_size*(cut_counts[d])
    end

    if l<length(grid.index)
        # This is true if we don't need to interpolate all dimensions because we're on a boundary:
        return grid.index[1:l]::Vector{Int}, grid.weight[1:l]::Vector{Float64}
    end
    grid.index::Vector{Int}, grid.weight::Vector{Float64}
end

function interpolants(grid::SimplexGrid, x::Vector)

    weight = grid.weight
    index  = grid.index

    x_p = grid.x_p # residuals
    ihi = grid.ihi # indicies of cuts above point
    ilo = grid.ilo # indicies of cuts below point
    n_ind = grid.n_ind

    cut_counts = grid.cut_counts
    cuts = grid.cuts

    cut_i = 1

    for i = 1:dimensions(grid)
        # find indicies of coords if match
        coord = x[i]
        lasti = cut_counts[i]+cut_i-1
        ii = cut_i
        # check bounds, snap to closest if out
        if coord <= cuts[ii]
            ihi[i] = ii
            ilo[i] = ii
            x_p[i] = 0.0
        elseif coord >= cuts[lasti]
            ihi[i] = lasti
            ilo[i] = lasti
            x_p[i] = 0.0
        else
            # increment through cut points if in bounds
            while cuts[ii] < coord
                ii += 1
            end
            # if on cut assign cut indecies
            if cuts[ii] == coord
                ilo[i] = ii
                ihi[i] = ii
                x_p[i] = 0.0
            else
                # if between cuts assign lo and high indecies and translate
                ilo[i] = ii-1
                ihi[i] = ii
                lo = cuts[ilo[i]]
                hi = cuts[ihi[i]]
                x_p[i] = (x[i] - lo) / (hi - lo)
            end
        end
        cut_i = cut_i + cut_counts[i]
    end

    # initialize sort indecies
    for i = 1:length(n_ind); n_ind[i] = i; end
    # sort translated and scaled x values
    sortperm!(n_ind, x_p, rev=true)
    x_p = x_p[n_ind]
    n_ind = n_ind - 1

    # get weight
    for i = 1:(length(x_p)+1)
        if i == 1
            weight[i] = 1 - x_p[i]
        elseif i == length(x_p)+1
            weight[i] = x_p[i-1]
        else
            weight[i] = x_p[i-1] - x_p[i]
        end
    end

    # get indecies
    fill!(index, 0)
    i_index = 0
    for i = 1:(length(x_p)+1)
        siz = 1
        ct = 0
        good_count = 1
        if i > 1
            i_index = i_index + 2^(n_ind[i-1])
        end
        for k = 1:length(x)
            onHi = ((i_index & good_count) > 0)
            good_count <<= 1
            if onHi
                index[i] += (ihi[k] - 1 - ct) * siz
            else
                index[i] += (ilo[k] - 1 - ct) * siz
            end
            siz = siz*cut_counts[k]
            ct += cut_counts[k]
        end
        index[i] += 1
    end

    weight = weight ./ sum(weight)

    return index::Vector{Int}, weight::Vector{Float64}
end

function interpolants(knn::KnnGrid, x::Vector, dist_metric::PreMetric = Euclidean())
    @assert length(x) == GridInterpolations.dimensions(knn)
    
    # Get list of possible k "nearest" values in each dimension
    possibilities_by_dim = [] 
    num_row = length(knn.cutPoints) # one row per dimension
    num_col = 1 # a column for every combination of k-values on each dimension
                # (but note that k-values can sometimes come back with length k+1) 
    for idx in 1:length(x)
        this_dimension = knn.cutPoints[idx] 
        p = get_nearest_k_in_dim(this_dimension, x[idx], knn.k)
        num_col *= length(p)
        push!(possibilities_by_dim, p)
    end

    # Convert the possibilities into a matrix that we can feed 
    # to a distance function
    possibilities = zeros(num_row, num_col)
    i = 1
    for p in product(possibilities_by_dim...) #############################################6 as killer of speed
        possibilities[:,i] = collect(p) 
        i += 1
    end

    # Get relative distances and their best order
    distances = colwise(dist_metric, possibilities, x)
    best_order = sortperm(distances)  ####################################### also killer of speed

    indices = zeros(Int, knn.k)
    for col in 1:knn.k
        ind_of_current_best = best_order[col]
        indices[col] = x2ind(knn,possibilities[:,ind_of_current_best])
    end
    weights = ones(knn.k)/knn.k
    
    return (indices::Vector{Int}, weights::Vector{Float64})
end


#################### sortperm! is included in Julia v0.4 ###################

using Base.Order # for sortperm!, should be availiable in v 0.4
using Base.Sort # for sortperm!
const DEFAULT_UNSTABLE = QuickSort

function sortperm!{I<:Integer}(x::Vector{I}, v::AbstractVector; alg::Algorithm=DEFAULT_UNSTABLE,
lt::Function=isless, by::Function=identity, rev::Bool=false, order::Ordering=Forward)
    sort!(x, alg, Perm(ord(lt,by,rev,order),v))
end

end # module




