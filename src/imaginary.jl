export imaginary_projection


"""
    imaginary_projection(f; alg=Greedy(), grid=..., options...)

Compute the imaginary_projection of `f` with the given algorithm `alg`
and a 2D or 3D grid. `alg` can be
* `Greedy()`
* `Simple()`
but you probably only want to use `Greedy()`.

## Examples

```julia
@polyvar x y

g = x^4 + im * x^3 - x^2 * y^2 + 3x^3 - 2im*x*y^2 + (4-2im)*x + 0.5*y^2 + 1.5
grid = Grid2D(xlims=(-4, 4), ylims=(-5, 5), res=(1200, 1200))
# This uses the `Greedy()` algorithm
imaginary_projection(g, grid=grid)
```

These algorithms are all approximations of the imaginary projection 𝐼(𝑓) based on a grid. This basically
applies the membership test for different grid points.
The grid can be passed explicitly, otherwise it will be computed based on a heuristic.

* `membership_options=[MembershipTestOptions()](@ref)`: The options for the membership test.
* `npasses=1`: The `Greedy()` algorithm can make multiple passes to improve the quality.
* `test_domain`: A tuple `(xmin, xmax, ymin, ymax, [zmin, zmax])` from which start values for the membership test are drawn.
"""
function imaginary_projection(f::MP.AbstractPolynomial; kwargs...)
    if MP.nvariables(f) == 2
        F = ImaginaryFiber2D(f)
    elseif MP.nvariables(f) == 3
        F = ImaginaryFiber3D(f)
    else
        throw(error("Currently only 2 and 3-dimensional PolynomialAmoebas are supported."))
    end
    imaginary_projection(F, f; kwargs...)
end

function imaginary_projection(
        F::AbstractImaginaryFiber,
        f::MP.AbstractPolynomial;
        alg::AbstractGridAlgorithm=Greedy(),
        kwargs...)
    imaginary_projection(F, f, alg; kwargs...)
end

function imaginary_projection(F::AbstractImaginaryFiber{T, N}, f::MP.AbstractPolynomial,
    alg::AbstractGridAlgorithm;
    grid=default_grid(F, f),
    callback=_do_nothing,
    membership_options=MembershipTestOptions(),
    npasses=1,
    test_domain=default_test_domain(F, f),
    ) where {T, N}
    if !(typeof(grid) <: AbstractGrid{N})
        return throw(error("The grid is not of dimension $N"))
    end
    if N == 3
        generator = DomainStartValueGenerator3D(test_domain)
    else
        generator = DomainStartValueGenerator2D(test_domain)
    end

    _imaginary_projection(alg, F, f, grid, generator, membership_options, callback, npasses)
end

@inline default_grid(F::AbstractImaginaryFiber{T, 2}, f::MP.AbstractPolynomial) where T = throw(error("You have to provide a `Grid2D` via grid=..."))
@inline default_grid(F::AbstractImaginaryFiber{T, 3}, f::MP.AbstractPolynomial) where T = throw(error("You have to provide a `Grid3D` via grid=..."))

function default_test_domain(F::AbstractImaginaryFiber{T, 2}, f::MP.AbstractPolynomial) where T
    try
        return amoeba_carcase_domain_heuristic(f, factor=3.0)
    catch err
        return (-20, 20, -20, 20)
    end
end
function default_test_domain(F::AbstractImaginaryFiber{T, 3}, f::MP.AbstractPolynomial) where T
    return (-20, 20, -20, 20, -20, 20)
end

function _imaginary_projection(::Simple, F, f::MP.AbstractPolynomial{<:Real}, grid, gen, options, callback, npasses)
    # We can exploit the fact that all solutions occur in conjugated pairs
    cache = falses(size(grid))
    avoid_check(B, k) = symmetry(B, k, grid, cache)

    simple_grid(F, grid, gen, avoid_check, options=options)
end

function symmetry(B, k, grid, cache)
    y = grid[k]
    sym_sub = grid[0 .- y]
    sym_k = _sub2ind(size(grid), sym_sub...)
    if 1 ≤ sym_k ≤ length(grid) && cache[sym_k]
        return true, B[sym_k]
    else
        cache[k] = true
        return false, false
    end
end

function _imaginary_projection(::Simple, F, f, grid, gen, options, callback, npasses)
    simple_grid(F, grid, gen, options=options)
end


function _imaginary_projection(::Greedy, F::ImaginaryFiber2D, f, grid, gen, options, callback, npasses)
    B, start_values = imaginary_by_amoeba(f, grid, options=options)
    initial_queue, queued = initial_neighbour_queue(B.data)
    B = greedy_grid_memorized(F, grid, gen, initial_queue, start_values, deepcopy(B), queued)
    for _ =2:npasses
        initial_queue, queued = initial_neighbour_queue(B.data)
        B = greedy_grid_memorized(F, grid, gen, initial_queue, start_values, deepcopy(B), queued)
    end

    B
end

function _imaginary_projection(::Greedy, F::ImaginaryFiber3D, f, grid, gen, options, callback, npasses)
    z = MP.variables(f)[end]
    B = empty(grid)
    start_values = Array{SVector{3, Float64}, 3}(undef, size(grid))
    grid_2d = Grid2D(grid.xrange, grid.yrange)
    for (k, z_k) in enumerate(grid.zrange)
        k % 2 == 0 && continue
        f_z_k = MP.subs(f, z=>z_k)
        B_z_k, start_values_z_k = imaginary_by_amoeba(f_z_k, grid_2d, options=options)
        B.data[:,:,k] = B_z_k.data
        start_values[:,:,k] = map(x -> SVector(x[1], x[2], 0.0), start_values_z_k)
    end
    initial_queue, queued = initial_neighbour_queue(B.data)
    B = greedy_grid_memorized(F, grid, gen, initial_queue, start_values, B, queued)

    B
end

function imaginary_by_amoeba(f, grid; options=MembershipTestOptions())
    start_values = Matrix{SVector{2, Float64}}(undef, size(grid)...)
    B = Bitmap2D(grid)
    m, n = size(B)
    function start_value_cb(result, k, bitmap)
        θ = result.solution
        w = bitmap.grid[k]
        r1, r2 = exp.(w)
        sθ = sin.(θ)
        cθ = cos.(θ)
        i, j = grid[(r1 * sθ[1], r2 * sθ[2])]
        if 1 ≤ i ≤ m && 1 ≤ j ≤ n
            B[i, j] = true
            start_values[i, j] = SVector(r1 * cθ[1], r2 * cθ[2])
        end
    end

    res = min(200, div(max(size(grid)...), 2))
    amoeba(f, alg=Greedy(), resolution=res, callback=start_value_cb)
    while !any(B.data)
        res *= 2
        amoeba(f, alg=Greedy(), resolution=res, callback=start_value_cb)
    end


    return B, start_values
end
