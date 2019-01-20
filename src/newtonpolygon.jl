export newtonpolygon, vertices, vertexindices, lattices, facets, subdivision

"""
    newtonpolygon(p::AbstractPolynomial; lowerhull=false)

Construct a [`NewtonPolygon`](@ref) from the support of the polynomial `p`. `lowerhull` indicates
whether the lower or upper convex hull should be used to construct the regular subdivision
of the Newton polygon.
"""
function newtonpolygon(p::MP.AbstractPolynomial; kwargs...)
    @assert MP.nvariables(p) == 2 "Expected a bivariate polynomial."
    lattices = map(t -> SVector{2}(MP.exponents(t)), MP.terms(p))
    coeffs = realcoefficients(p)
    return newtonpolygon(lattices, coeffs; kwargs...)
end
function newtonpolygon(p::MP.AbstractPolynomial{<:Complex})
    @assert MP.nvariables(p) == 2 "Expected a bivariate polynomial."
    lattices = map(t -> SVector{2}(MP.exponents(t)), MP.terms(p))
    return newtonpolygon(lattices)
end

function newtonpolygon(lattices::Vector{SVector{2, Int}})
    vertices, facets = vertices_facets(lattices)
    NewtonPolygon(lattices, vertices, facets, Vector{Vector{Int}}())
end
function newtonpolygon(lattices::Vector{SVector{2, Int}}, coefficients::AbstractVector{<:Real}; lowerhull=false)
    vertices, facets = vertices_facets(lattices)
    try
        subdivision = newtonsubdivision(lattices, coefficients, lowerhull)
        return NewtonPolygon(lattices, vertices, facets, subdivision)
    catch err
        subdivision = [vertices]
        return NewtonPolygon(lattices, vertices, facets, subdivision)
    end

end

function newtonpolygon(exponents::Matrix)
    @assert size(exponents, 1) == 2 "Expected a bivariate polynomial"
    lattices = [SVector(exponents[1, j], exponents[2, j]) for j=1:size(exponents, 2)]
    newtonpolygon(lattices)
end

function vertices_facets(lattices)
    vertices = Vector{Int}()
    try
        hull = convexhull([lattices[i][j] for i=1:length(lattices), j=1:2])
        vertices = hull.vertices
    catch err
        throw(AssertionError("Newton polygon is not fulldimensional"))
        # we have a degenerate case...
        if length(lattices) == 2
            vertices = [1, 2]
        elseif length(lattices) == 1
            vertices = [1]
        elseif length(lattices) > 2
            # the points are coplanar
            vertices = [1, length(lattices)]
        end
    end

    # vertices are in counter clockwise orders
    facets = Vector{SVector{2,Int}}()
    n = length(vertices)
    for i = 1:n
        push!(facets, SVector(vertices[i], vertices[i%n + 1]))
    end

    vertices, facets
end

"""
    newtonsubdivision(lattices, coefficients::AbstractVector{<:Real})

compute the subdivion of the newton polygon. This can throw if the lattices
are not generic enough
"""
function newtonsubdivision(lattices, coeffs::AbstractVector{<:Real}, lowerhull=false)
    nterms = length(lattices)
    subdivision = Vector{Vector{Int}}()
    if nterms == 3
        push!(subdivision, [1, 2, 3])
    elseif nterms > 3
        liftedvertices = Matrix{Float64}(undef, nterms, 3)
        for i in eachindex(lattices)
            liftedvertices[i, 1] = lattices[i][1]
            liftedvertices[i, 2] = lattices[i][2]
            liftedvertices[i, 3] = coeffs[i]
        end


        hull = convexhull(liftedvertices)
        simplices = Vector{Vector{Int}}()
        equations = Vector{Vector{Float64}}()
        for (i, simplex) in enumerate(hull.simplices)
            if lowerhull ?  hull.equations[i,3] < 0.0 : hull.equations[i,3] > 0.0
                push!(simplices, simplex)
                push!(equations, hull.equations[i,:])
            end
        end

        # Problem is now that the scipy always returns a triangulation
        # Therefore we have to manually put the simplices together again
        n = length(simplices)

        to_get_ordered = Vector{Vector{Int}}()
        handled = falses(n)
        for i=1:n
            if handled[i]
                continue
            end
            eq = equations[i]
            s = simplices[i]
            merged_something = false
            for j=i+1:n
                if handled[j]
                    continue
                end
                if equations[j] == eq
                    append!(s, simplices[j])
                    merged_something = true
                    handled[j] = true
                end
            end
            if merged_something
                push!(to_get_ordered, unique(s))
            else
                push!(subdivision, s)
            end
        end

        for x in to_get_ordered
            h = convexhull([lattices[i][j] for i in x, j=1:2])
            push!(subdivision, x[h.vertices])
        end
    else
        error("Newton polygon is not full dimensional")
    end
    subdivision
end

"""
    lattices(newton_polygon)

The lattice points of the Newton polygon
"""
lattices(np::NewtonPolygon) = np.lattices

"""
    vertexindices(newton_polygon)

The indices of the vertices of the Newton polygon.
"""
vertexindices(np::NewtonPolygon) = np.vertices


"""
    vertices(newton_polygon)

The coordinates of the vertices of the Newton polygon.
"""
vertices(np::NewtonPolygon) = np.lattices[np.vertices]


"""
    facets(newton_polygon)

The facets of the Newton polygon. Each facet is a `SVector{2}` with the coordinates
of the start and end vertex of the facet.
"""
facets(np::NewtonPolygon) = map(f -> SVector(np.lattices[f[1]], np.lattices[f[2]]), np.facets)

iscontained(np::NewtonPolygon, p) = iscontained(vertices(np), p)
on_boundary(np::NewtonPolygon, p) = on_boundary(vertices(np), p)

"""
    subdivision(newton_polygon)

The regular subdivision of the Newton polygon. Each simplex is a vector with the coordinates
of the vertices of the simplex.
"""
function subdivision(np::NewtonPolygon)
    map(np.subdivision) do s
        np.lattices[s]
    end
end

@recipe function plot(newt::NewtonPolygon; showsubdivision=false)
    xs, ys = unzip(vertices(newt))
    xlims --> (minimum(xs) - 1, maximum(xs) + 1)
    ylims --> (minimum(ys) - 1, maximum(ys) + 1)
    border --> nothing
    legend --> nothing
    aspect_ratio --> :equal
    @series begin
        seriestype := :shape
        fill --> :transparent
        linewidth --> 1
        [xs; first(xs)], [ys; first(ys)]
    end

    if showsubdivision
        for simplex in subdivision(newt)
            @series begin
                seriestype := :shape
                fill --> :transparent
                linewidth --> 1
                unzip([simplex..., simplex[1]])
            end
        end
    end

    @series begin
        seriestype := :scatter
        markercolor --> :black
        markerstrokewidth --> 0
        markersize --> 6
        unzip(lattices(newt))
    end
end
