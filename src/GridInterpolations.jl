module GridInterpolations

using StaticArrays
using Printf
using LinearAlgebra

export AbstractGrid, RectangleGrid, SimplexGrid, dimensions, length, label, ind2x, ind2x!, interpolate, maskedInterpolate, interpolants, vertices

abstract type AbstractGrid{D} end # D is the dimension

mutable struct RectangleGrid{D} <: AbstractGrid{D}
    cutPoints::Vector{Vector{Float64}}
    cut_counts::Vector{Int}
    cuts::Vector{Float64}

    function RectangleGrid{D}(cutPoints...) where {D}
        cut_counts = Int[length(cutPoints[i]) for i = 1:length(cutPoints)]
        cuts = vcat(cutPoints...)
        myCutPoints = Array{Vector{Float64}}(undef, length(cutPoints))
        numDims = length(cutPoints)
        @assert numDims == D
        for i = 1:numDims
            if length(Set(cutPoints[i])) != length(cutPoints[i])
                error(@sprintf("Duplicates cutpoints are not allowed (duplicates observed in dimension %d)", i))
            end
            if !issorted(cutPoints[i])
                error("Cut points must be sorted")
            end
            myCutPoints[i] = cutPoints[i]
        end
        return new(myCutPoints, cut_counts, cuts)
    end
end

RectangleGrid(cutPoints...) = RectangleGrid{length(cutPoints)}(cutPoints...)

mutable struct SimplexGrid{D} <: AbstractGrid{D}
    cutPoints::Vector{Vector{Float64}}
    cut_counts::Vector{Int}
    cuts::Vector{Float64}
    x_p::Vector{Float64} # residuals
    ihi::Vector{Int} # indices of cuts above point
    ilo::Vector{Int} # indices of cuts below point
    n_ind::Vector{Int}

    function SimplexGrid{D}(cutPoints...) where {D}
        cut_counts = Int[length(cutPoints[i]) for i = 1:length(cutPoints)]
        cuts = vcat(cutPoints...)
        myCutPoints = Array{Vector{Float64}}(undef, length(cutPoints))
        numDims = length(cutPoints)
        @assert numDims == D
        for i = 1:numDims
            if length(Set(cutPoints[i])) != length(cutPoints[i])
                error(@sprintf("Duplicates cutpoints are not allowed (duplicates observed in dimension %d)", i))
            end
            if !issorted(cutPoints[i])
                error("Cut points must be sorted")
            end
            myCutPoints[i] = cutPoints[i]
        end
        x_p = zeros(numDims) # residuals
        ihi = zeros(Int, numDims) # indicies of cuts above point
        ilo = zeros(Int, numDims) # indicies of cuts below point
        n_ind = zeros(Int, numDims)
        return new(myCutPoints, cut_counts, cuts, x_p, ihi, ilo, n_ind)
    end
end

SimplexGrid(cutPoints...) = SimplexGrid{length(cutPoints)}(cutPoints...)

Base.length(grid::RectangleGrid) = prod(grid.cut_counts)
Base.size(grid::RectangleGrid) = Tuple(grid.cut_counts)
Base.length(grid::SimplexGrid) = prod(grid.cut_counts)

dimensions(grid::AbstractGrid{D}) where {D} = D
Base.ndims(grid::AbstractGrid{D}) where {D} = D

label(grid::RectangleGrid) = "multilinear interpolation grid"
label(grid::SimplexGrid) = "simplex interpolation grid"

# showall returns a valid constructor incantation - it will be called when repr is called on a grid
function Base.show(io::IO, grid::AbstractGrid)
    if get(io, :compact, false)
        print(io, "$(typeof(grid)) with $(length(grid)) points")
    else
        print(io, "$(typeof(grid))(")
        for v in grid.cutPoints
            show(io, v)
            print(io, ',')
        end
        print(io, ')')
    end
end

function ind2x(grid::AbstractGrid, ind::Int)
    ndims = dimensions(grid)
    x = Array{Float64}(undef, ndims)
    ind2x!(grid, ind, x)
    x::Array{Float64}
end

function ind2x!(grid::AbstractGrid, ind::Int, x::AbstractArray)
    # Populates x with the value at ind.
    # In-place version of ind2x.
    # Example:
    #   rgrid = RectangleGrid([2,5],[20,50])
    #   x = [0,0]
    #   ind2x!(rgrid,4,x)  # x now contains [5,50]
    #   @show x            # displays [5,50]
    ndims = dimensions(grid)
    stride = grid.cut_counts[1]
    for i = 2:ndims-1
        stride *= grid.cut_counts[i]
    end

    for i = (ndims-1):-1:1
        rest = rem(ind - 1, stride) + 1
        x[i+1] = grid.cutPoints[i+1][div(ind - rest, stride)+1]
        ind = rest
        stride = div(stride, grid.cut_counts[i])
    end
    x[1] = grid.cutPoints[1][ind]
    nothing
end


# masked interpolation ignores points that are masked
function maskedInterpolate(grid::AbstractGrid, data::DenseArray, x::AbstractVector, mask::BitArray{1})
    index, weight = interpolants(grid, x)
    val = zero(eltype(data))
    totalWeight = zero(eltype(data))
    for i = 1:length(index)
        if mask[index[i]]
            continue
        end
        val += data[index[i]] * weight[i]
        totalWeight += weight[i]
    end
    return val / totalWeight
end

interpolate(grid::AbstractGrid, data::Matrix, x::AbstractVector) = interpolate(grid, map(eltype(data), data[:]), x)

function interpolate(grid::AbstractGrid, data::DenseArray, x::AbstractVector)
    index, weight = interpolants(grid, x)
    v = zero(eltype(data))
    for (i, data_ind) in enumerate(index)
        v += data[data_ind] * weight[i]
    end
    return v
end

function interpolants(grid::RectangleGrid, x::AbstractVector)
    if any(isnan, x)
        throw(DomainError("Input contains NaN!"))
    end
    cut_counts = grid.cut_counts
    cuts = grid.cuts

    # Reset the values in index and weight:
    num_points = 2^dimensions(grid)
    index = MVector{num_points, Int}(undef)
    index2 = MVector{num_points, Int}(undef)
    weight = MVector{num_points, eltype(x)}(undef)
    weight2 = MVector{num_points, eltype(x)}(undef)

    # Note: these values are set explicitly because we have not verified that the logic below is independent of the initial values. See discussion in PR #47. These can be removed if it can be proved that the logic is independent of the initial values.
    index .= 1
    index2 .= 1
    weight .= zero(eltype(weight))
    weight2 .= zero(eltype(weight2))

    weight[1] = one(eltype(weight))
    weight2[1] = one(eltype(weight2))

    l = 1
    subblock_size = 1
    cut_i = 1
    n = 1
    for d = 1:length(x)
        coord = x[d]
        lasti = cut_counts[d] + cut_i - 1
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
                i_lo, i_hi = (ii - 1), ii
            end
        end

        # the @inbounds are needed below to prevent allocation
        if i_lo == i_hi
            for i = 1:l
                @inbounds index[i] += (i_lo - cut_i) * subblock_size
            end
        else
            low = (1 - (coord - cuts[i_lo]) / (cuts[i_hi] - cuts[i_lo]))
            for i = 1:l
                @inbounds index2[i] = index[i] + (i_lo - cut_i) * subblock_size
                @inbounds index2[i+l] = index[i] + (i_hi - cut_i) * subblock_size
            end
            @inbounds index[:] = index2
            for i = 1:l
                @inbounds weight2[i] = weight[i] * low
                @inbounds weight2[i+l] = weight[i] * (1 - low)
            end
            @inbounds weight[:] = weight2
            l = l * 2
            n = n * 2
        end
        cut_i = cut_i + cut_counts[d]
        subblock_size = subblock_size * (cut_counts[d])
    end

    v = min(l, length(index))
    return view(SVector(index), 1:v), view(SVector(weight), 1:v)
end

function interpolants(grid::SimplexGrid, x::AbstractVector)

    weight = MVector{dimensions(grid) + 1,eltype(x)}(undef)
    index = MVector{dimensions(grid) + 1,Int}(undef)

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
        lasti = cut_counts[i] + cut_i - 1
        ii = cut_i
        # check bounds, snap to closest if out
        if coord <= cuts[ii]
            ihi[i] = ii
            ilo[i] = ii
            x_p[i] = zero(eltype(x))
        elseif coord >= cuts[lasti]
            ihi[i] = lasti
            ilo[i] = lasti
            x_p[i] = zero(eltype(x))
        else
            # increment through cut points if in bounds
            while cuts[ii] < coord
                ii += 1
            end
            # if on cut assign cut indecies
            if cuts[ii] == coord
                ilo[i] = ii
                ihi[i] = ii
                x_p[i] = zero(eltype(x))
            else
                # if between cuts assign lo and high indecies and translate
                ilo[i] = ii - 1
                ihi[i] = ii
                lo = cuts[ilo[i]]
                hi = cuts[ihi[i]]
                x_p[i] = (x[i] - lo) / (hi - lo)
            end
        end
        cut_i = cut_i + cut_counts[i]
    end

    # initialize sort indecies
    for i = 1:length(n_ind)
        n_ind[i] = i
    end
    # sort translated and scaled x values
    sortperm!(n_ind, x_p, rev=true) ############################################# killer of speed
    x_p = x_p[n_ind]
    n_ind = n_ind .- 1

    # get weight
    for i = 1:(length(x_p)+1)
        if i == 1
            weight[i] = 1 - x_p[i]
        elseif i == length(x_p) + 1
            weight[i] = x_p[i-1]
        else
            weight[i] = x_p[i-1] - x_p[i]
        end
    end

    # get indices
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
            siz = siz * cut_counts[k]
            ct += cut_counts[k]
        end
        index[i] += 1
    end

    weight = weight ./ sum(weight)

    return SVector(index), SVector(weight)
end

"Return a vector of SVectors where the ith vector represents the vertex corresponding to the ith index of grid data."
function vertices(grid::AbstractGrid)
    n_dims = dimensions(grid)
    mem = Array{Float64,2}(undef, n_dims, length(grid))

    for idx = 1:length(grid)
        this_idx::Int = idx - 1

        # Get the correct index into each dimension
        # and populate vertex index with corresponding cut point
        for j = 1:n_dims
            cut_idx::Int = this_idx % grid.cut_counts[j]
            this_idx = div(this_idx, grid.cut_counts[j])
            mem[j, idx] = grid.cutPoints[j][cut_idx+1]
        end
    end

    #=
    This relies on the memory layout of Matrix to stay the same, so is a
    possible source of future errors. However, it is documented
    (http://juliaarrays.github.io/StaticArrays.jl/stable/pages/
    api.html#Arrays-of-static-arrays-1), and tests should catch these errors.
    =#
    return reshape(reinterpret(SVector{n_dims,Float64}, mem), (length(grid),))
end



#################### sortperm! is included in Julia v0.4 ###################

using Base.Order # for sortperm!, should be availiable in v 0.4
using Base.Sort # for sortperm!
const DEFAULT_UNSTABLE = QuickSort

function sortperm!(x::Vector{I}, v::AbstractVector; alg::Algorithm=DEFAULT_UNSTABLE,
    lt::Function=isless, by::Function=identity, rev::Bool=false, order::Ordering=Forward) where {I<:Integer}
    sort!(x, alg, Perm(ord(lt, by, rev, order), v))
end

function Base.iterate(iter::RectangleGrid, state::Int64=1)
    return state <= length(iter) ? (ind2x(iter, state), state + 1) : nothing
end

Base.getindex(grid::RectangleGrid, key::CartesianIndex) = ind2x(grid, LinearIndices(Dims((grid.cut_counts...,)))[key])
Base.getindex(grid::RectangleGrid, indices...) = ind2x(grid, LinearIndices(Dims((grid.cut_counts...,)))[indices...])

end # module
