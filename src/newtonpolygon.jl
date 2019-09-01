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

function newtonpolygon(lattice_points::Vector{SVector{2, Int}}, lift::Union{Nothing,AbstractVector{<:Real}}=nothing; lowerhull=false)
    P = Polyhedra.polyhedron(Polyhedra.vrep(lattice_points), CDDLib.Library())
    Polyhedra.removevredundancy!(P)
    if Polyhedra.dim(P) != 2
        throw(AssertionError("Newton polygon is not fulldimensional"))
    end

    vpts = sort_counterclockwise(collect(Polyhedra.points(P)))
    vertices = map(vpts) do v
        for (i,p) in enumerate(lattice_points)
            p == v && return i
        end
    end

    # facets
    idx = Int[]
    for p in vpts
        for vid in vertices
            if p == lattice_points[vid]
                push!(idx, vid)
                break
            end
        end
    end
    # iterate reverse to obtain counter-clockwise order
    facets = [SVector(idx[i], idx[i+1]) for i in 1:length(idx)-1]
    push!(facets, SVector(idx[end], idx[1]))

    # subdivison
    subdivision = Vector{Int}[]
    if lift !== nothing
        lifted_pts = map((p, w) -> [p;w], lattice_points, lift)
        lifted_P = Polyhedra.polyhedron(Polyhedra.vrep(lifted_pts), CDDLib.Library())
        Polyhedra.removevredundancy!(lifted_P)
        if Polyhedra.dim(lifted_P) == 2
            push!(subdivision, copy(vertices))
        else
            hspaces = Polyhedra.halfspaces(lifted_P)
            for (idx, H) in zip(eachindex(hspaces), hspaces)
                if (lowerhull ? H.a[end] < 0 : H.a[end] > 0)
                    ipts = Polyhedra.incidentpoints(lifted_P, idx)
                    pts = sort_counterclockwise(map(p -> SVector(p[1], p[2]), ipts))
                    push!(subdivision, map(pts) do p
                        for (i, v) in enumerate(lattice_points)
                            p[1] == v[1] && p[2] == v[2] && return i
                        end
                    end)
                end
            end
        end
    end

    return NewtonPolygon(lattice_points, vertices, facets, subdivision)
end

function getsemihull(ps::Vector{PT}, sign_sense, counterclockwise, yray = nothing) where PT
    hull = PT[]
    if length(ps) == 0
        return hull
    end
    prev = sign_sense == 1 ? first(ps) : last(ps)
    cur = prev
    for j in (sign_sense == 1 ? (2:length(ps)) : ((length(ps)-1):-1:1))
        while prev != cur && counterclockwise(cur - prev, ps[j] - prev) >= 0
            cur = prev
            pop!(hull)
            if !isempty(hull)
                prev = last(hull)
            end
        end
        if yray !== nothing && counterclockwise(ps[j] - cur, yray) >= 0
            break
        else
            push!(hull, cur)
            prev = cur
            cur = ps[j]
        end
    end
    push!(hull, cur)
    hull
end

function sort_counterclockwise(ps)
    sort!(ps, by = first)
    counterclockwise(p1, p2) = dot(cross([p1; 0], [p2; 0]), [0, 0, 1])
    top = getsemihull(ps,  1, counterclockwise)
    bot = getsemihull(ps, -1, counterclockwise)
    if !isempty(top) && !isempty(bot)
        @assert top[end] == bot[1]
        pop!(top)
    end
    if !isempty(bot) && !isempty(top)
        @assert bot[end] == top[1]
        pop!(bot)
    end
    sort!(top, rev=true)
    sort!(bot)
    [top; bot]
end

function newtonpolygon(exponents::Matrix)
    @assert size(exponents, 1) == 2 "Expected a bivariate polynomial"
    lattices = [SVector(exponents[1, j], exponents[2, j]) for j=1:size(exponents, 2)]
    newtonpolygon(lattices)
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

"""
    convexhull_vertices(line::Vector{Tuple{Int, Int}})

Returns a list of indices of the vertices of the convex hull.
"""
function convexhull_vertices(points::Vector{NTuple{2, T}}) where {T<:Real}
    convexhull_vertices(SVector.(points))
end
function convexhull_vertices(points::Vector{SVector{2, T}}) where {T<:Real}
    P = Polyhedra.polyhedron(Polyhedra.vrep(points), CDDLib.Library())
    Polyhedra.removevredundancy!(P)
    if Polyhedra.dim(P) != 2
        throw(ArgumentError("points are not full dimensional"))
    end
    vpts = sort_counterclockwise(collect(Polyhedra.points(P)))
    vertices = map(vpts) do v
        for (i,p) in enumerate(points)
            p == v && return i
        end
    end
    vertices
end
