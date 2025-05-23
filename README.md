# GridInterpolations

[![Build Status](https://github.com/sisl/GridInterpolations.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/sisl/GridInterpolations.jl/actions/workflows/CI.yml)
[![Coverage Status](https://coveralls.io/repos/sisl/GridInterpolations.jl/badge.svg)](https://coveralls.io/r/sisl/GridInterpolations.jl)

This package performs multivariate interpolation on a rectilinear grid. At the moment, it provides implementations of multilinear and simplex interpolation. As of benchmarks in December 2016, multilinear interpolation performs fastest and with the most accuracy.

The following image visualizes grid-based interpolation in two dimensions, with shape of interpolater for (−0.3,0.8) inscribed. The small dots reﬂect the interpolation's estimate for sin(x)+2cos(y)+sin(5xy), which is the underlying reward function approximated by the large dot lattice.

![Illustration of performance of multilinear and simplex interpolation methods](sampleInterpolation.png)

For a description of multilinear and simplex interpolation see: Scott Davies, _Multidimensional Triangulation and Interpolation for Reinforcement Learning_, Advances in Neural Information Processing Systems, Cambridge, MA: MIT Press, 1997. [pdf](http://papers.nips.cc/paper/1229-multidimensional-triangulation-and-interpolation-for-reinforcement-learning.pdf)

There are some related packages, such as [Interpolations.jl](https://github.com/tlycken/Interpolations.jl).  

## Installation

Start Julia and run the following command:

```julia
Pkg.add("GridInterpolations")
```

## Usage

To use the GridInterpolations module, begin your code with

```julia
using GridInterpolations
```

## Interpolation

Create two-dimensional interpolation grids, a data array, and a point of interest:
```julia
grid = RectangleGrid([0., 0.5, 1.],[0., 0.5, 1.])  	# rectangular grid
sGrid = SimplexGrid([0., 0.5, 1.],[0., 0.5, 1.])	# simplex grid
nGrid = NearestGrid([0., 0.5, 1.], [0., 0.5, 1.])     # nearest-neighbor grid 
gridData = [8., 1., 6., 3., 5., 7., 4., 9., 2.]   	# vector of value data at each cut
x = [0.25, 0.75]  									# point at which to perform interpolation
```

Perform interpolation on the rectangular grid:
```julia
julia> interpolate(grid,gridData,x)
5.25
```

Or interpolate on the simplex grid:
```julia
julia> interpolate(sGrid,gridData,x)
6.0
```

Or interpolate on the nearest-neighbor grid:
```julia
julia> interpolate(nGrid, gridData, x)
9.0
```

Compute interpolants for the grids:
```julia
julia> sGrid = SimplexGrid([0., 0.5, 1.],[0., 0.5, 1.])
[[0.0,0.5,1.0],[0.0,0.5,1.0]]

julia> interpolants(sGrid, x)
([4,5,8],[0.5,0.0,0.5])
```

Convert an index to a Grid coordinate:
```julia
julia> ind2x(grid, 3)
2-element Array{Float64,1}:
 1.0
 0.0
```

Number of vertices in the grid:
```julia
julia> length(grid)
9
```

Number of dimensions:
```julia
julia> dimensions(grid)
2
```

Multi-dimensional indexing using Cartesian coordinates:
```julia
julia> [grid[c] for c in CartesianIndices((3,3))]
3×3 Array{Array{Float64,1},2}:
 [0.0, 0.0]  [0.0, 0.5]  [0.0, 1.0]
 [0.5, 0.0]  [0.5, 0.5]  [0.5, 1.0]
 [1.0, 0.0]  [1.0, 0.5]  [1.0, 1.0]
```
or multi-dimensional indices
```julia
julia> grid[2,2]
2-element Array{Float64,1}:
 0.5
 0.5
```

Sequential iteration over grid points:
```julia
julia> for x in grid
           # do stuff
       end
```

### Auto differentiation

Autodiff packages are also supported:
```julia
grid_data = [8.0, 1.0, 6.0, 3.0, 5.0, 7.0, 4.0, 9.0, 2.0]
x = [0.25, 0.75]
grid = RectangleGrid([0.0, 0.5, 1.0], [0.0, 0.5, 1.0])

using ForwardDiff
f(x::Vector) = interpolate(grid, grid_data, x)
ForwardDiff.gradient(f, x)
```

## Limitations

`RectangleGrid` will error in high-dimensional domains (above about 15). In these cases `SimplexGrid` should be used.

## Credits

Contributors to this package include Maxim Egorov, Eric Mueller, and Mykel Kochenderfer
