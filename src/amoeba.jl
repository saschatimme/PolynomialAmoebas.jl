export amoeba

"""
    amoeba(f; alg=Polygonal(), options...)

Compute the amoeba of `f` with the given algorithm `alg` which can be
* `Polygonal()`
* `Greedy()`
* `ArchTrop()`
* `Simple()`

but you probably only want to use `Polygonal()` or `Greedy()`.

Depending on the given algorithm the function takes different keyword arguments.

## Examples

```julia
@polyvar x y

# This uses the `Polygonal()` algorithm
amoeba(x^2 + y^2 + 1)

# We only want a crude approximation
amoeba(x^2 + y^2 + 1, accuracy=0.1)

# Use the `Greedy()` algorithm.
amoeba(x^2 + y^2 + 1, alg=Greedy())

# Use the `Greedy()` algorithm with a custom grid.
grid = Grid2D(xlims=(-5, 5), ylims=(-4, 4), res=(500, 400))
amoeba(x^2 + y^2 + 1, alg=Greedy(), grid=grid)

# Use the `Greedy()` algorithm with the default domain but a higher resolution
amoeba(x^2 + y^2 + 1, alg=Greedy(), resolution=800)
```

## `Polygonal()`
This algorithm computes an approximation of the amoeba ``\\mathcal{A}(f)`` in the provided `domain` Ω
by computing a set of polygons 𝒫 with |𝒫| ⊂ Ω such that the union of these polygons approximates ``\\mathcal{A}(f)`` ∩ Ω
from the outside.

The possible (optional) arguments are
* `domain`: A tuple in the form `(xmin, xmax, ymin, ymax)` which defines a section Ω
for which the amoeba ``\\mathcal{A}(f)`` is computed. This domain has to be such that the intersection
Ω ∩ ``\\mathcal{A}(f)`` still captures the correct topology of ``\\mathcal{A}(f)``.
* `accuracy=0.01`: The maximal allowed error |𝒫 - ``\\mathcal{A}(f)`` ∩ Ω|. Note that we only compute an upper limit of the error.
The algorithms stops if the given accuracy is reached.
* `spine`: This algorithm needs the spine of ``\\mathcal{A}(f)``.
* `minimal_component_size=0.01`: The minimal size of the components of the complement. This is only used if no spine
is passed explicitly.
* `iterations=2000`: The maximal number of iterations.
* `vertices_accuracy=accuracy*1e-3`: During the algorithm we approximate points on the boundary of ``\\mathcal{A}(f)``. This
is the accuracy with which we compute them. Note hat this influences the minimal error of the approximation
and it should always be some magnitudes smaller than `accuracy`.
* `membership_options=[MembershipTestOptions()](@ref)`: As a subroutine a membership test is used.


## `Greedy()`, `Simple()`, `ArchTrop()`

These algorithms are all approximations of the amoeba ``\\mathcal{A}(f)`` based on a grid. This basically
applies the membership test for different grid points. `Greedy()` is the fastest and
`Simple()` the slowest.
The grid can be passed explicitly, otherwise it will be computed based on a heuristic.

* `resolution=600`: The resolution of the grid if not passed explicitly.
* `grid`: If passed explicitly this grid is used, otherwise it will be computed based on a heuristic.
* `membership_options=[MembershipTestOptions()](@ref)`: The options for the membership test.
"""
function amoeba end

function amoeba(f::MP.AbstractPolynomial; kwargs...)
    if MP.nvariables(f) == 2
        F = AmoebaFiber2D(f)
    elseif MP.nvariables(f) == 3
        F = AmoebaFiber3D(f)
    else
        throw(error("Currently only 2 and 3-dimensional PolynomialAmoebas are supported."))
    end
    amoeba(F, f; kwargs...)
end

function amoeba(
        F::AbstractAmoebaFiber,
        f::MP.AbstractPolynomial;
        alg::AbstractAlgorithm=default_algorithm(F),
        kwargs...)
    amoeba(F, f, alg; kwargs...)
end

default_algorithm(F::AbstractAmoebaFiber{T, 2}) where T = Polygonal()
default_algorithm(F::AbstractAmoebaFiber{T, 3}) where T = Greedy()

function amoeba(F::AmoebaFiber3D, f, alg::Polygonal; kwargs...)
    throw(error("The `Polygonal` algorithm is only supported for bivariate polynomials."))
end
function amoeba(
        F::AmoebaFiber2D{T},
        f::MP.AbstractPolynomial,
        alg::Polygonal;
        domain=default_amoeba_domain(f),
        accuracy=1e-2,
        vertices_accuracy = accuracy * 1e-3,
        iterations = 2000,
        membership_options = MembershipTestOptions(),
        minimal_component_size = 0.01,
        spine = Spine2D(f, minimal_component_size=minimal_component_size, domain=domain),
        report_progress=false
        ) where {T}

    C = Covering(spine, domain...)
    M = CoveringFitter(C, F;
        desired_maximal_error=accuracy, vertices_accuracy=vertices_accuracy,
        maximal_iterations=iterations, options=membership_options)
    fit!(M, report_progress=report_progress)

    polygons = cartesian_polygons(M.covering)
    error_bound = M.approximated_global_error

    PolygonalAmoeba(spine, polygons, float.(domain), error_bound, M.iteration)
end

function amoeba(F::AbstractAmoebaFiber{T, N}, f::MP.AbstractPolynomial,
    alg::AbstractGridAlgorithm;
    resolution=600,
    grid=default_grid(F, f, resolution),
    callback=_do_nothing,
    membership_options=MembershipTestOptions()
    ) where {T, N}
    if !(typeof(grid) <: AbstractGrid{N})
        return throw(error("The grid is not of dimension $N"))
    end

    _amoeba(alg, F, f, grid, TorusStartValueGenerator{N}(), membership_options, callback)
end

@inline default_grid(F::AbstractAmoebaFiber{T, 2}, f::MP.AbstractPolynomial, r) where T = Grid2D(f, r)
@inline default_grid(F::AbstractAmoebaFiber{T, 3}, f::MP.AbstractPolynomial, r) where T = throw(error("You have to provide a `Grid3D` via grid=..."))


function default_amoeba_domain(f::MP.AbstractPolynomial; factor=1.5, aspect_ratio=:default)
    try
        return amoeba_carcase_domain_heuristic(f, factor=factor, aspect_ratio=aspect_ratio)
    catch err
        throw(error("You have to pass an explicit domain."))

    end
end

_amoeba(::Simple, F, f, grid, gen, options, callback) = simple_grid(F, grid, gen, options=options)
function _amoeba(::ArchTrop, F::AmoebaFiber2D, f, grid, gen, options, callback)
    nb = archimedean_neighbourhood(f, grid).data
    avoid_check = (B, k) -> (!nb[k], false)
    simple_grid(F, grid, gen, avoid_check, options=options)
end

function _amoeba(::Greedy, F::AmoebaFiber2D, f, grid, gen, options, callback)
    initial_queue = findall(grid_contour(f, grid).data)

    greedy_grid(F, grid, gen, initial_queue, options=options, callback=callback)
end

function _amoeba(::Greedy, F::AmoebaFiber3D, f, grid, gen, options, callback)
    n1, n2, n3 = size(grid)

    z = MP.variables(f)[end]
    k0 = div(n3, 3)
    f_z0 = MP.subs(f, z => exp(grid.zrange[k0]))
    grid_k0 = Grid2D(grid.xrange, grid.yrange)
    contour_f_z0 = grid_contour(f_z0, grid_k0).data
    initial_queue = map(findall(contour_f_z0)) do ij
        CartesianIndex(Tuple(ij)..., k0)
    end

    greedy_grid(F, grid, gen, initial_queue, options=options, callback=callback)
end
