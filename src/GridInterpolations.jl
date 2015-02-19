module GridInterpolations

export AbstractGrid, RectangleGrid, SimplexGrid, dimensions, length, show, ind2x, ind2x!, interpolate, maskedInterpolate, interpolants

abstract AbstractGrid

type RectangleGrid <: AbstractGrid
    cutPoints::Vector{Vector{Float64}}
    cut_counts::Vector{Int}
    cuts::Vector{Float64}
    index::Vector{Int}
    weight::Vector{Float64}
    index2::Vector{Int}
    weight2::Vector{Float64}
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
end

Base.length(grid::RectangleGrid) = prod(grid.cut_counts)

Base.length(grid::SimplexGrid) = prod(grid.cut_counts)

dimensions(grid::RectangleGrid) = length(grid.cut_counts)

dimensions(grid::SimplexGrid) = length(grid.cut_counts)

function RectangleGrid(cutPoints::Vector{Float64}...)
    cut_counts = Int[length(cutPoints[i]) for i = 1:length(cutPoints)]
    cuts = vcat(cutPoints...)
    myCutPoints = Array(Vector{Float64}, length(cutPoints))
    for i = 1:length(cutPoints)
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
    RectangleGrid(myCutPoints, cut_counts, cuts, index, weight, index2, weight2)
end

function SimplexGrid(cutPoints::Vector{Float64}...)
    cut_counts = Int[length(cutPoints[i]) for i = 1:length(cutPoints)]
    cuts = vcat(cutPoints...)
    myCutPoints = Array(Vector{Float64}, length(cutPoints))
    for i = 1:length(cutPoints)
        myCutPoints[i] = cutPoints[i]
    end
    numDims = length(cutPoints)
    index = zeros(Int, numDims+1) # d+1 points for simplex
    weight = zeros(Float64, numDims+1)
    x_p = zeros(numDims) # residuals
    ihi = zeros(Int, numDims) # indicies of cuts above point
    ilo = zeros(Int, numDims) # indicies of cuts below point
    n_ind = zeros(Int, numDims)
    SimplexGrid(myCutPoints, cut_counts, cuts, index, weight, x_p, ihi, ilo, n_ind)
end

Base.show(io::IO, grid::AbstractGrid) = show(io, grid.cutPoints)

function ind2x(grid::AbstractGrid, ind::Int)
    ndims = dimensions(grid)
    x = Array(Float64, ndims)
    ind2x!(grid, ind, x)
    x::Array{Float64}
end

function ind2x!(grid::AbstractGrid, ind::Int, x::Array{Float64})
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

# masked interpolation ignores points that are masked
function maskedInterpolate(grid::AbstractGrid, data::Vector{Float64}, x::Vector{Float64}, mask::BitArray{1})
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

interpolate(grid::AbstractGrid, data::Array, x::Vector{Float64}) = interpolate(grid, float64(data[:]), x)

function interpolate(grid::AbstractGrid, data::Vector{Float64}, x::Vector{Float64})
    index, weight = interpolants(grid, x)
    dot(data[index], weight)
end

function interpolants(grid::RectangleGrid, x::Vector{Float64})
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

function interpolants(grid::SimplexGrid, x::Vector{Float64})

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

#################### sortperm! is included in Julia v0.4 ###################

using Base.Order # for sortperm!, should be availiable in v 0.4
using Base.Sort # for sortperm!
const DEFAULT_UNSTABLE = QuickSort

function sortperm!{I<:Integer}(x::Vector{I}, v::AbstractVector; alg::Algorithm=DEFAULT_UNSTABLE,
lt::Function=isless, by::Function=identity, rev::Bool=false, order::Ordering=Forward)
    sort!(x, alg, Perm(ord(lt,by,rev,order),v))
end

end # module
