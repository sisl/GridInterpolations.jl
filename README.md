# GridInterpolations
## Usage:

To use the GridInterpolations module, begin your code with

```julia
using GridInterpolations
```

### Interpolation

Create a rectangular and a simplex interpolation grids in two dimensions, a data array, and a point of interest:
```julia
grid = RectangleGrid([0., 0.5, 1.],[0., 0.5, 1.])  	# rectangular grid
sGrid = SimplexGrid([0., 0.5, 1.],[0., 0.5, 1.])	# simplex grid
gridData = [8. 1. 6. 3. 5. 7. 4. 9. 2.]            	# value data at each cut
x = [0.25, 0.75]  									# point at which to perform interpolation
```

Perform the interpolation on the rectangular grid:
```julia
julia> interpolate(grid,gridData,x)
5.25
```

Or interpolate on the simplex grid:
```julia
julia> interpolate(sGrid,gridData,x)
6.0
```

Both simplexMagicTest20.txt and rectangularMagicTest20.txt are needed in order to successfully execute runtests.jl.

## Credits

The main authors of this package are Mykel Kochenderfer, Eric Mueller, and Maxim Egorov.