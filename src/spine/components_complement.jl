export order, point

"""
    order(cc::ComponentComplement)

The order of the given component of the complement.
"""
order(cc::ComponentComplement) = cc.order

"""
    point(cc::ComponentComplement)

A point in the given component of the complement.
"""
point(cc::ComponentComplement) = cc.x

"""
    isbounded(cc::ComponentComplement)

Returns whether the component of the complement `cc` is bounded.
"""
isbounded(cc::ComponentComplement) = cc.bounded

function Base.show(io::IO, cc::ComponentComplement{N}) where N
    print(io,
        "ComponentComplement{$N} with order ",
        tuple(cc.order...),
        cc.bounded ? " (bounded)" : " (unbounded)",)
end

@recipe function plot(CC::Vector{ComponentComplement}; show_orders=true, show_markers=!show_orders)
    legend --> false
    markerstrokewidth --> 0.0
    if show_markers
        markersize --> 6.0
    else
        markersize --> 0.0
    end
    markercolor --> :orange
    @series begin
        st := :scatter
        points = point.(CC)
        xs, ys = unzip(points)

        if show_orders
            series_annotations := map(order.(CC)) do ord
                "($(ord[1]),$(ord[2]))"
            end
        end

        xs, ys
    end
end


"""
    components_complement(amoeba_approximation, f; nsamples=512)

Compute a all components of the complement of the amoeba of `f` from the given
`amoeba_approximation` `Bitmap2D`.
"""
function components_complement(M::Bitmap2D, f::MP.AbstractPolynomial{T}; kwargs...) where T
    components_complement(M, AmoebaFiber2D(f), SP.Polynomial(promote_type(T, Float64), f); kwargs...)
end
function components_complement(M::Bitmap2D, F::AmoebaFiber2D, f::SP.Polynomial{T}; nsamples=1024) where {T}
    newt = newtonpolygon(SP.exponents(f))
    boundaries = complement_boundaries(M, f, newt, nsamples)
    components = Vector{ComponentComplement{2}}()
    # we avoid some allocations in `order` with this
    working_vec = zeros(complex(T), 2)
    gen = TorusStartValueGenerator{2}()
    options = MembershipTestOptions(maxiters=200, ntries=50)
    for (α, boundary) in boundaries
        try
            c = point_in_complement(boundary, α, f, M.grid)
            bounded = !on_boundary(newt, α)
            # We validate that the points is definitely not in the amoeba
            if membershiptest(F, (c[1], c[2]), gen, options).successfull
                continue
            end
            push!(components, ComponentComplement{2}(α, c, bounded))
        catch err
        end
    end

    components, boundaries
end

function point_in_complement(component_boundary, component_order, f::SP.Polynomial{T}, grid) where T
    working_vec = zeros(complex(T), 2)
    k = 0
    while k < 10
        try
            chull = convexhull_vertices(component_boundary)
            c = poly_centroid(component_boundary[chull], grid)
            ord = round.(Int, real.(order(f, c, working_vec, 1024)))
            if ord == component_order
                return c
            else
                c = boundary_average(component_boundary, grid)
                ord = round.(Int, real.(order(f, c, working_vec, 1024)))
                if ord == component_order
                    return c
                else
                    # We discard all vertices of the convex hull as outliers
                    # (One of them is an outlier)
                    deleteat!(component_boundary, chull)
                end
            end
        catch err
            # An error indicates that our complement component is not 2-dimensional
            # but rather a line. I.e. it is sufficient to consider the
            # average of all (unique) grid points on the line
            points = sort!(unique(component_boundary))
            if length(points) == 1
                c = grid[points[1]]
            else
                m = div(length(points), 2)
                c = grid[points[m]...]
            end
            ord = round.(Int, real.(order(f, c, working_vec, 1024)))
            if ord == component_order
                return c
            else
                throw(ErrorException("Cannot locate point in component with order $component_order."))
            end
        end
    end
    throw(ErrorException("Cannot locate point in component with order $component_order."))
end

"""
    boundary_average(boundary, grid)

Compute the average of the contour
"""
function boundary_average(boundary::Vector{NTuple{2, Int}}, grid::Grid2D)
    mean(p -> grid[p], boundary)
end

function centroid(b::Vector{NTuple{2, Int}}, grid::Grid2D)
    try
        return poly_centroid(b[convexhull_vertices(b)], grid)
    catch
        # An error indicates that our complement component is not 2-dimensional
        # but rather a line. I.e. it is sufficient to consider the
        # average of all (unique) grid points on the line
        points = sort!(unique(b))
        if length(points) == 1
            return grid[points[1]]
        else
            m = div(length(points), 2)
            return grid[points[m]]
        end
    end
end


function poly_centroid(boundary::Vector{NTuple{2, Int}}, grid::Grid2D)
    A = 0.0
    cx = 0.0
    cy = 0.0

    n = length(boundary)
    i = n
    j = 1

    while j ≤ n
        xi, yi = grid[boundary[i]]
        xi1, yi1 = grid[boundary[j]]
        a = xi * yi1 - xi1 * yi
        cx += (xi + xi1) * a
        cy += (yi + yi1) * a
        A += a
        i = j
        j += 1
    end
    # @show A, length(boundary)
    (1 / (3A)) .* SVector(cx, cy)
end


"""
    complement_boundaries(M::Bitmap2D, f, newt)

Compute the boundaries of all complement components (as grid indices).
Note that this is noisy. There can be duplicate components usually arising
from small additional components (with a boundary of length 5).
"""
function complement_boundaries(M::Bitmap2D, f::SP.Polynomial, newt, nsamples)
    # we will add a padding of 1 around the image to get closed contours only
    # Therefore we have to add two points to the x and y-coordinates
    # x
    boundaries = contour_boundaries(M)
    sort!(boundaries, lt=((a, b) -> length(a) < length(b)), rev=true)
    split_with_order_map(boundaries, M.grid, f, newt, nsamples)
end

function contour_boundaries(M)
    # we will add a padding of 1 around the image to get closed contours only
    # Therefore we have to add two points to the x and y-coordinates
    # x
    xcoords = xcoordinates(M.grid)
    xspace = [xcoords[1] - step(xcoords)]
    append!(xspace, xcoords)
    push!(xspace, xcoords[end] + step(xcoords))

    ycoords = ycoordinates(M.grid)
    yspace = [ycoords[1] - step(ycoords)]
    append!(yspace, ycoords)
    push!(yspace, ycoords[end] + step(ycoords))

    values = floatpadd(M.data)
    lines = Contour.contour(xspace, yspace, values, 1e-30) |> Contour.lines
    boundaries = map(lines) do line
        xs, ys = Contour.coordinates(line)
        line_grid_points = Vector{NTuple{2, Int}}()
        for (x, y) in zip(xs, ys)
            i, j = gridpoint(M.grid, x, y)
            i = clamp(i, 1, size(M.grid, 1))
            j = clamp(j, 1, size(M.grid, 2))
            if !M[i,j]
                push!(line_grid_points, (i, j))
            end
        end
        line_grid_points
    end
end

"""
    split_with_order_map(boundaries, grid, f, newt)

Split the boundaries by using the order map.
"""
function split_with_order_map(boundaries, G, f::SP.Polynomial{T}, newt, nsamples) where T
    splitted_boundaries = Vector{eltype(boundaries)}()
    splitted_orders = Vector{SVector{2, Int}}()
    # we need to handle the "jump" between different orders
    index_bookkeeping = Dict{SVector{2, Int}, Int}()
    last_order = (-1, -1)
    last_index = 0
    working_vec = zeros(complex(T), 2)
    for boundary in boundaries
        for p in boundary
            imag_ord = order(f, G[p...], working_vec, nsamples)
            ord = round.(Int, real.(imag_ord))
            if norm(ord - imag_ord, Inf) > 1e-4
                continue
            end
            if ord == last_order
                index = last_index
            elseif haskey(index_bookkeeping, ord)
                index = index_bookkeeping[ord]
            elseif iscontained(newt, ord)
                push!(splitted_boundaries, Vector{Tuple{Int, Int}}())
                push!(splitted_orders, ord)
                index = length(splitted_boundaries)
                push!(index_bookkeeping, ord => index)
            else
                continue
            end

            push!(splitted_boundaries[index], p)

            last_index = index
            last_order = ord
        end
    end
    zip(splitted_orders, splitted_boundaries)
end

"""
    order(f, w)

Compute the order ``(v_1, v_2)`` of the given complements.

The ``j``-th order of ``f`` at ``w=(w_1,w_2)`` can be computed by the number of zeros
(minus the number of poles) of the one-variable Laurent-polynomial
```math
    u ↦ f(c_1^u^{δ_{1j}}, c_2^u^{δ_{2j}})
```
inside the unit circle ``|u|=1``, where ``(c1, c2)`` is any vector with ``Log(c1, c2) = w`` [^1].

In order to compute this we can use the
[argument principle](https://en.m.wikipedia.org/wiki/Argument_principle).
For this we have to rewrite the integral (here only for the first component of the order)
```math
\frac{1}{2πi}∫_{|u|=1}\frac{f'(e^{w_1}u, e^{w_2})}{f(e^{w_1}u, e^{w_2})}du
```
With a change of variables we get
```math
\frac{1}{2πi}∫_0^{2π}\frac{f'(e^{w_1}e^{iθ}, e^{w_2})}{f(e^{w_1}e^{iθ}, e^{w_2})}ie^{iθ}dθ
```
[^1]: Forsberg, Mikael, Mikael Passare, and August Tsikh. "Laurent determinants and arrangements of hyperplane PolynomialAmoebas." Advances in mathematics 151.1 (2000): 45-70.
"""
function order(f::SP.Polynomial, w, working_vec::AbstractVector, nsamples=128)
    v1 = 0.0im
    v2 = 0.0im
    h = im * 2π / nsamples
    c1, c2 = exp.(w)
    for k=0:nsamples-1
        u = exp(h * k)
        # g1(θ) = f(c1 * exp(iθ), c2)
        # g1'(θ) = c1 * f_x(c1 * exp(iθ), c2)
        # p1(θ) = g1'(θ) * i * exp(iθ) / g1(θ)
        # v1 += p1(θ)
        c1u = c1 * u
        x = SVector(c1u, c2)
        f_x = SP.gradient(f, x)[1]
        val = SP.evaluate(f, x)
        v1 += c1u * f_x / val

        # g2(θ) = f(c1, c2 * exp(iθ))
        # g2'(θ) = c2 * f_y(c1, c2 * exp(iθ))
        # p2(θ) = g2'(θ) * i * exp(iθ) / g2(θ)
        # v2 += p1(θ)
        c2u = c2 * u
        x = SVector(c1, c2u)
        f_y = SP.gradient(f, x)[2]
        val = SP.evaluate(f, x)
        v2 += c2u * f_y / val
    end
    v1 = v1 / (nsamples)
    v2 = v2 / (nsamples)
    SVector{2}(v1, v2)
end





"Add a rectangle of `true`s around the `BitArray`"
function floatpadd(A::BitArray)
    m, n = size(A)
    M = ones(m + 2, n + 2)
    for i=1:m, j=1:n
        M[i+1, j+1] = float(A[i, j])
    end
    M
end
