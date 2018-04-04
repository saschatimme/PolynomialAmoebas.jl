export contour


"""
    contour(f; options...)

Compute the contour 𝐶(𝑓) of the amoeba ``\\mathcal{A}(f)``.

## Example
```julia
@polyvar x y

contour(x^2+y^2+1)

# custom domain
contour(x^2+y^2+1, domain=(-5, 5, -5, 5))
```

## Options
* `domain`: A tuple in the form `(xmin, xmax, ymin, ymax)` which defines a section Ω
for which the contour 𝐶(𝑓) is computed.
* `membership_options`: The options for the membership test
* `res=(600,600)`: The resolution with which starting points are sampled
* `samples_off_axis=2*MP.maxdegree(p)^2`: The number of sample points per off-axis.
"""
function contour(f::MP.AbstractPolynomial;
    domain=amoeba_carcase_domain_heuristic(f, factor=2.0),
    res::NTuple{2, Int}=(600, 600),
    samples_off_axis::Int=2*MP.maxdegree(f),
    membership_options=MembershipTestOptions(ntries=10, maxiters=30, tol=1e-8))

    grid = Grid2D(domain..., res...)

    _contour2D(ContourFiber2D(f), grid, samples_off_axis, membership_options)
end

# only internally used
"""
    grid_contour(F::ContourFiber2D, grid::Grid2D; kwargs...)

Instead of returning a list of contour points it fills a bitmap with the nearest contour
points.
"""
function grid_contour(p::MP.AbstractPolynomial; res=600, kwargs...)
    grid_contour(p, Grid2D(p, res); kwargs...)
end

function grid_contour(p::MP.AbstractPolynomial, grid::Grid2D;
    every_main_axis::Int=10, samples_off_axis::Int=8, options=MembershipTestOptions(ntries=3, maxiters=40, tol=1e-8))

    fiber = ContourFiber2D(p)

    _contour2D_grid(fiber, grid, every_main_axis, samples_off_axis, options)
end
