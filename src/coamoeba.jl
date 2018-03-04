export coamoeba

"""
    coamoeba(f; alg=Greedy(), options...)

Compute the coamoeba of `f` with the given algorithm `alg` which can be
* `Greedy()`
* `Coarse()`
* `Simple()`
but you probably only want to use `Greedy()` or `Coarse`.

The coamoeba is embedded in ``[0, 2π)^n`` where ``n`` is either 2 or 3 depending
on the polynomial. The algorithm approximate the coamoeba on a grid representing this
embedding.

## Example
```julia
@polyvar x y

# This uses the `Greedy()` algorithm
coamoeba(x^2 + y^2 + 1)
```

## Optional arguments
* `resolution=600`: The resolution of the grid if not passed explicitly.
* `membership_options=[MembershipTestOptions()](@ref)`: The options for the membership test.
* `test_domain`: Tuple `(xmin, xmax, ymin, ymax)` resp. `(xmin, xmax, ymin, ymax, zmin, zmax)`
from which start values for the membership test are drawn if necessary.
"""
function coamoeba(f::MP.AbstractPolynomial; alg::AbstractAlgorithm=Greedy(), kwargs...)
    if MP.nvariables(f) == 2
        fiber = CoamoebaFiber2D(f)
    elseif MP.nvariables(f) == 3
        fiber = CoamoebaFiber3D(f)
    else
        throw(error("Polynomial has $(MP.nvariables(f)) variables. Expected 2 or 3."))
    end
    coamoeba(fiber, f, alg; kwargs...)
end


function coamoeba(F::CoamoebaFiber2D, f::MP.AbstractPolynomial, alg;
    test_domain=default_coamoeba_domain(f),
    resolution=200,
    membership_options=MembershipTestOptions(ntries=20))

    grid = Grid2D((0.0, 2.0π), resolution)
    generator = DomainStartValueGenerator2D(test_domain...)

    _coamoeba(alg, F, f, grid, generator, membership_options)
end

function coamoeba(F::CoamoebaFiber3D, f::MP.AbstractPolynomial, alg;
    test_domain=(-5, 5, -5, 5, -5, 5),
    resolution=21,
    membership_options=MembershipTestOptions(ntries=20))

    grid = Grid3D((0.0, 2.0π), resolution)
    generator = DomainStartValueGenerator3D(test_domain...)

    _coamoeba(alg, F, f, grid, generator, membership_options)
end


function default_coamoeba_domain(f::MP.AbstractPolynomial)
    try
        return amoeba_carcase_domain_heuristic(f, factor=3.0)
    catch err
        return (-5, 5, -5, 5)
        # throw(error("You have to pass an explicit domain."))
    end
end


function _coamoeba(::Coarse, fiber::CoamoebaFiber2D, f, grid, generator, options)
    coamoeba = _coarse_coamoeba(f, grid)
    Coamoeba(coamoeba)
end

function _coamoeba(::Simple, fiber::CoamoebaFiber2D, f, grid, generator, options)
    coarse = Bitmap2D(_coarse_coamoeba(f, grid), grid)
    Coamoeba(simple_grid(fiber, grid, generator, _always_check, coarse, options=options).data)
end

function _coamoeba(::Greedy, fiber::CoamoebaFiber2D, f, grid, generator, options)
    start_values = fill(SVector(NaN, NaN), size(grid)...)
    B = Bitmap2D(_coarse_coamoeba(f, grid), grid)
    initial_queue, queued =initial_neighbour_queue(B.data)

    result = greedy_grid_memorized(fiber, grid, generator, initial_queue, start_values, B, queued)

    Coamoeba{2}(result.data)
end

function _coamoeba(::Simple, fiber, f, grid, generator, options)
    Coamoeba(simple_grid(fiber, grid, generator, options=options).data)
end

function _coamoeba(::Coarse, F::CoamoebaFiber3D, f, grid, gen, options)
    n1, n2, n3 = size(grid)
    coamoeba = falses(n1, n2, n3)
    x, y, z = MP.variables(f)
    for k0 = 1:n3
        f_z0 = MP.subs(f, z => cis(grid.zrange[k0]))
        grid_k0 = Grid2D(grid.xrange, grid.yrange)
        coarse_f_z0 = _coarse_coamoeba(f_z0, grid_k0)

        coamoeba[:,:,k0] .= coarse_f_z0
    end

    for j0 = 1:n2
        f_y0 = MP.subs(f, y => cis(grid.yrange[j0]))
        grid_j0 = Grid2D(grid.xrange, grid.zrange)
        coarse_f_y0 = _coarse_coamoeba(f_y0, grid_j0, @view coamoeba[:,j0,:])

        coamoeba[:,j0,:] .= coarse_f_y0
    end
    for i0 = 1:n1
        f_x0 = MP.subs(f, x => cis(grid.xrange[i0]))
        grid_i0 = Grid2D(grid.yrange, grid.zrange)
        coarse_f_x0 = _coarse_coamoeba(f_x0, grid_i0, @view coamoeba[i0,:,:])

        coamoeba[i0,:,:] .= coarse_f_x0
    end

    Coamoeba(coamoeba)
end

function _coamoeba(::Greedy, fiber::CoamoebaFiber3D, f, grid, generator, options)
    coarse = _coamoeba(Coarse(), fiber, f, grid, generator, options).data
    queued = copy(coarse)
    initial_queue = Int[]
    for k = 1:length(coarse)
        if coarse[k]
            addneighbours!(initial_queue, queued, k)
        end
    end
    bitmap = Bitmap3D(coarse, grid)
    Coamoeba(greedy_grid(fiber, grid, generator, initial_queue, bitmap, queued, options=options).data)
end

function _coarse_coamoeba(f, grid, coamoeba = falses(size(grid)))
    cvs = coefficients_vertices(f, newtonpolygon(f))
    partition = partition_work(length(grid))
    Threads.@threads for tid = 1:Threads.nthreads()
        r = partition[tid]
        for k in r
            if !coamoeba[k]
                coamoeba[k] = winding_number(cvs, grid[k]) != 0
            end
        end
    end
    coamoeba
end

"""
    coefficients_vertices(f, newt::NewtonPolygon)

Return a list of tuples `(c_k/|c_k|, v_k)` where `v_k` is the k-th vertex
of the Newton polygon and `c_k` the corresponding coefficient. The vertices
are in counter-clockwise orientation.
"""
function coefficients_vertices(f, newt::NewtonPolygon)
    map(vertices(newt)) do v
        for t in MP.terms(f)
            if MP.exponents(t) == v
                c = MP.coefficient(t)
                return (c / abs(c), v)
            end
        end
    end
end

const _inv_2π = inv(2π)

function winding_number(cvs, θ)
    # k_1 = length(cvs)
    c, vert = cvs[end]
    @inbounds prev_v = c * cis(muladd(vert[1], θ[1], vert[2] * θ[2]))
    w = 0.0im
    for k = 1:length(cvs)
        c, vert = cvs[k]
        @inbounds v = c * cis(muladd(vert[1], θ[1], vert[2] * θ[2]))
        w -= log(prev_v / v)
        prev_v = v
    end
    w *= _inv_2π
    round(Int, w.im)
end


function Base.show(io::IO, A::Coamoeba{2})
    n, m = size(A.data)
    println(io, "Coamoeba of size ", n, "×", m)
end
function Base.show(io::IO, A::Coamoeba{3})
    n1, n2, n3 = size(A.data)
    println(io, "Coamoeba of size ", n1, "×", n2, "×", n3)
end

function Base.show(io::IO, mime::MIME"text/html", A::Coamoeba)
    if inIJulia() && plotsdefined()
        show(io, mime, Main.Plots.plot(A))
    else
        show(io, A)
    end
end

@recipe function plot(M::Coamoeba{2})
    data = [ M.data M.data
             M.data M.data]
    G = Grid2D(-2π, 2π, -2π, 2π, size(data, 1), size(data, 2))
    @series begin
        xticks := ([-2π,-π, 0, π, 2π], ["-2\\pi", "-\\pi", "0", "\\pi", "2\\pi"])
        yticks := ([-2π,-π, 0, π, 2π], ["-2\\pi", "-\\pi", "0", "\\pi", "2\\pi"])
        Bitmap2D(data, G)
    end
end

@recipe function plot(M::Coamoeba{3})
    grid = Grid3D(xlims=(0,2π), ylims=(0,2π), zlims=(0,2π), res=size(M.data))
    @series begin
        xticks := ([-2π,-π, 0, π, 2π], ["-2\\pi", "-\\pi", "0", "\\pi", "2\\pi"])
        yticks := ([-2π,-π, 0, π, 2π], ["-2\\pi", "-\\pi", "0", "\\pi", "2\\pi"])
        zticks := ([-2π,-π, 0, π, 2π], ["-2\\pi", "-\\pi", "0", "\\pi", "2\\pi"])
        Bitmap3D(M.data, grid)
    end
end
