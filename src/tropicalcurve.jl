export TropicalCurve, vertices, segments, halfrays, archtrop_neighbourhoods_masked


"""
    TropicalCurve(tropicalpoly)
    TropicalCurve(tropicalpoly, newton_polygon::NewtonPolygon)

Computes a Tropical planar curve. A curce is represented as a vector of `vertices`, a vector
of `segments` between vertices and `halfrays` starting from a node in a specify direction.

You can visualize a `TropicalCurve` `curve` by using `Plots.jl` via `plot(curve)`.


    TropicalCurve(realpoly)

The polynomial will first be converted to a tropical one via [tropicalpolynomial](@ref)
and then the curve will be computed as usual.
"""
function TropicalCurve(polynomial::AbstractTropicalPolynomial, newt)
    coeffs = realcoefficients(polynomial)
    vertices = _vertices(newt, coeffs)
    segments = _segments(newt)
    halfrays = _halfrays(newt)

    TropicalCurve(polynomial, vertices, segments, halfrays)
end
TropicalCurve(p::AbstractTropicalPolynomial{<:Real}) = TropicalCurve(p, newtonpolygon(p))
TropicalCurve(p::MP.AbstractPolynomial{<:Real}) = TropicalCurve(p, tropicalpolynomial(p))

function Base.show(io::IO, curve::TropicalCurve)
    println(io, "TropicalCurve from $(curve.polynomial)")
    println(io, "vertices: ", curve.vertices)
    println(io, "segments: ", curve.segments)
    println(io, "halfrays: ", curve.halfrays)
end


function _vertices(newt::NewtonPolygon, coeffs)
    vertices = map(newt.subdivision) do simplex
        # even if our subdivision has more than 3 vertices,
        # the first 3 are sufficient to compute the node
        k, l, m = simplex
        # by duality the minimum should be attained at least three times, i.e.
        #   max(c_{i_l, j_l} + i_l x + j_l y) = max(c_{i_k, j_k} + i_k x + j_k y)
        #   max(c_{i_m, j_m} + i_m x + j_m y) = max(c_{i_k, j_k} + i_k x + j_k y)
        i_k, j_k = newt.lattices[k]
        c_k = coeffs[k]
        i_l, j_l = newt.lattices[l]
        c_l = coeffs[l]
        i_m, j_m = newt.lattices[m]
        c_m = coeffs[m]

        A = SMatrix{2,2, Float64}(i_l - i_k, i_m - i_k, j_l - j_k, j_m - j_k)
        b = SVector(c_k - c_l, c_k - c_m)
        node = inv(A) * b
    end
    vertices
end

function _segments(newt::NewtonPolygon)
    n = length(newt.subdivision)
    segments = Vector{SVector{2, Int}}()
    for i = 1:n-1, j = i+1:n
        if commonedge(newt.subdivision[i], newt.subdivision[j])
            push!(segments, SVector(i,j))
        end
    end
    segments
end

function _halfrays(newt::NewtonPolygon)
    rays = Vector{Tuple{Int64, SVector{2, Float64}}}()
    for (i, simplex) in enumerate(newt.subdivision)
        n = length(simplex)
        for k in 1:n
            a, b = simplex[k], simplex[k%n + 1]
            # we take the midpoint of a facet of the simplex and check whether it is on the
            # convex hull of the newton polygon
            classification, val = classifypoint(newt, 0.5 * newt.lattices[a] + 0.5 * newt.lattices[b])
            if classification == :hull
                push!(rays, (i, val))
            end
        end
    end
    rays
end

"""
    classifypoint(newt, p)

Check whether a points is on the convex hull of the Newton polygon, and if so return the
outer normal.
"""
function classifypoint(newt::NewtonPolygon, p::SVector{2, Float64})
    for facet in newt.facets
        i, j = facet
        a = newt.lattices[i]
        b = newt.lattices[j]
        if distance(a, b, p) < 1e-6
            # we have a point on the convex hull
            # the facets are sorted counterclockwise,
            # i.e. to get an outer normal we have to rotate 90 degree clockwise
            return (:hull, normalize(SVector(b[2] - a[2], a[1] - b[1])))
        end
    end
    (:interior, p)
end


function commonedge(a::Vector{Int}, b::Vector{Int})
    # a = sort(a)
    # b = sort(b)

    m = length(a)
    n = length(b)
    for i in 1:m
        x = a[i]
        y = a[i%m+1]
        for j = 1:n
            if (b[j] == x && b[j%n+1] == y) || (b[j%n+1] == x && b[j] == y)
                return true
            end
        end
    end
    return false
end

"""
    vertices(tropicalcurve)

The vertices of the tropical curve.
"""
vertices(c::TropicalCurve) = c.vertices

"""
    segments(tropicalcurve)

The segments of the tropical curve. Each segment is a `SVector{2}` with start and end
coordinates.
"""
segments(c::TropicalCurve) = map(s -> c.vertices[s], c.segments)


"""
    halfrays(tropicalcurve)

The halfrays of the tropical curve. Each halfray is a tuple `(coordinate, direction)`.
"""
halfrays(c::TropicalCurve) = map(pair -> (c.vertices[pair[1]], pair[2]), c.halfrays)

"""
    curve2grid(f, grid, tropicalcurve)

Map the tropical curve `tropicalcurve` on the `Grid` `grid` and call function `f` for each grid
point  on that curve.
"""
function curve2grid(f::Function, grid::Grid2D, curve::TropicalCurve)
    for segment in curve.segments
        a, b = curve.vertices[segment]
        line(f, grid, nearestpixel(grid, a), nearestpixel(grid, b))
    end

    for (raynode, raydirection) in curve.halfrays
        a = curve.vertices[raynode]
        b = a + maxlength(a, limits(grid)...) * raydirection
        line(f, grid, nearestpixel(grid, a), nearestpixel(grid, b))
    end
    nothing
end


@recipe function plot(curve::TropicalCurve)
    @assert !isempty(curve.vertices) "Curve has no vertices"

    xmin, xmax, ymin, ymax = amoeba_carcase_domain_heuristic(curve)
    if !haskey(plotattributes, :xlims)
        xlims --> (xmin, xmax)
    else
        xmin, xmax = float.(plotattributes[:xlims])
    end
    if !haskey(plotattributes, :ylims)
        ylims --> (ymin, ymax)
    else
        ymin, ymax = float.(plotattributes[:ylims])
    end

    # grid --> false
    legend --> nothing
    aspect_ratio --> :equal
    linecolor --> Colors.colorant"#CD6155"
    linewidth --> 1.0



    rays = map(curve.halfrays) do x
        raynode, raydirection = x
        a = curve.vertices[raynode]
        b, _, _ = intersection_halfray_limits(a, normalize(raydirection), xmin, xmax, ymin, ymax)
        (a, b)
    end

    for segment in curve.segments
        a, b = curve.vertices[segment]
        @series begin
            aspect_ratio --> :equal
            fill := :transparent
            seriestype := :shape
            label := ""
            [a[1], b[1]], [a[2], b[2]]
        end
    end

    for (a, b) in rays
        @series begin
            aspect_ratio --> :equal
            fill := :transparent
            seriestype := :shape
            label := ""
            [a[1], b[1]], [a[2], b[2]]
        end
    end

end
