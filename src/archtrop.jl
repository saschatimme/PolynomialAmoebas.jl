export amoeba_carcase_domain_heuristic, archimedean_tropical_polynomial, archimedean_tropical_curve, archimedean_neighbourhood, archimedean_tropical_hypersurface

"""
    archimedean_tropical_polynomial(f)

Compute the Archimedean tropical polynomial of `f`.
"""
function archimedean_tropical_polynomial(f::MP.AbstractPolynomial{<:Union{Real, Complex}})
    mapcoefficients(c -> Tropical(logabs.(c)), f)
end

"""
    archimedean_tropical_curve(f)

Compute the curve associated to the Archimedean tropical polynomial of `f`.
"""
function archimedean_tropical_curve(f::MP.AbstractPolynomial{<:Union{Real, Complex}})
    MP.nvariables(f) == 2 || throw(error("Only bivariate polynomials are currently supported"))
    TropicalCurve(archimedean_tropical_polynomial(f))
end
archimedean_tropical_curve(f::AbstractTropicalPolynomial) = TropicalCurve(f)

"""
    archimedean_tropical_curve(f)

Compute the hypersurface associated to the Archimedean tropical polynomial of `f`.
"""
function archimedean_tropical_hypersurface(f::MP.AbstractPolynomial{<:Union{Real, Complex}})
    MP.nvariables(f) == 2 || throw(error("Only bivariate polynomials are currently supported"))
    archimedean_tropical_curve(f)
end


"""
    amoeba_carcase_domain_heuristic(f; factor=1.5, aspect_ratio=:default)

Compute the boundary `(xmin, xmax, ymin, ymax)` of domain Ω such that
Ω ∩ Amoeba(f) is probably the carcase of Amoeba(f).
"""
function amoeba_carcase_domain_heuristic(f::MP.AbstractPolynomial; kwargs...)
    amoeba_carcase_domain_heuristic(archimedean_tropical_curve(f); kwargs...)
end
function amoeba_carcase_domain_heuristic(curve::TropicalCurve; factor=1.5, aspect_ratio=:default)
    t = MP.nterms(curve.polynomial)
    r = factor * log(t - 1)

    xmin, xmax = Inf, -Inf
    ymin, ymax = Inf, -Inf
    for node in curve.vertices
        x, y = node
        xmin = min(xmin, x - r)
        xmax = max(xmax, x + r)
        ymin = min(ymin, y - r)
        ymax = max(ymax, y + r)
    end

    if aspect_ratio == :equal
        return equal_aspect_ratio(xmin, xmax, ymin, ymax)
    end

    xmin, xmax, ymin, ymax
end

"""
    archimedean_neighbourhood(archtrop, grid)

Compute the neighbourhood of `archtrop` in which the amoeba of `archtrop.polynomial` is
guaranteed to lie. See Avendaño, Kogan, Nisse, Rojas [^1] for details.

[^1]: Avendano, Martin, et al. "Metric Estimates and membership complexity
for archimedean amoebae and tropical hypersurfaces." arXiv preprint arXiv:1307.3681 (2013).
"""
function archimedean_neighbourhood(curve::TropicalCurve, grid)
    t = MP.nterms(curve.polynomial)
    r = log(t - 1)
    nb = archimedean_neighbourhood_rectangles(curve, r, limits(grid)...)
    M = falses(size(grid)...)
    for (px, py) in nb
        fillrect!(M, xgridpoint.(grid, px), ygridpoint.(grid, py))
    end
    Bitmap2D(M', grid)
end
function archimedean_neighbourhood(f::MP.AbstractPolynomial, grid)
    archimedean_neighbourhood(archimedean_tropical_curve(f), grid)
end

"Compute the neighbourhoods with radius `r` around the curve"
function archimedean_neighbourhood_rectangles(curve::TropicalCurve, r, xmin, xmax, ymin, ymax)
    # This doesn't really work if we the ray starts from outside the visible area
    nbs = Vector{NTuple{2, SVector{4,Float64}}}()
    for (raynode, raydirection) in curve.halfrays
        a = curve.vertices[raynode]
        b = a + maxlength(a, xmin, xmax, ymin, ymax) * raydirection
        push!(nbs, _rectline(a, b, r))
    end

    for segment in curve.segments
        a, b = curve.vertices[segment]
        # in degenerate cases it can happen that vertices coincide
        if norm(a-b) > 1e-8
            push!(nbs, _rectline(a, b, r))
        end
    end
    nbs
end
@inline function maxlength(a::SVector{2}, xmin, xmax, ymin, ymax)
    max(
        norm(a - SVector(xmin, ymin)),
        norm(a - SVector(xmin, ymax)),
        norm(a - SVector(xmax, ymin)),
        norm(a - SVector(xmax, ymax))
        )
end

"A rectangle centered around the line from `a` to `b` with height `2r`"
function _rectline(a, b, r)
    v1, v2 = normalize(a - b)
    d1 = -r*v2
    d2 = r*v1
    (SVector(a[1] + d1, b[1] + d1, b[1] - d1, a[1] - d1),
     SVector(a[2] + d2, b[2] + d2, b[2] - d2, a[2] - d2))
end

# function Base.show(io::IO, mime::MIME"text/html", A::ArchTrop)
#     if inIJulia() && plotsdefined()
#         show(io, mime, Main.Plots.plot(curve(A)))
#     else
#         show(io, A)
#     end
# end
# function Base.show(io::IO, A::ArchTrop)
#     #println(typeof(A), ":")
#     println(io, "ArchTrop of $(A.polynomial) is $(A.curve.polynomial)")
# end
