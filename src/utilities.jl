
_isone(x::T) where T = x == one(T)
#
function _ind2sub(a, i)
     i2s = CartesianIndices(a)
     Tuple(i2s[i])
end

function _sub2ind(a,i...)
      s2i = LinearIndices(a)
      s2i[i...]
end

function Base.isless(a::SVector{2}, b::SVector{2})
    if a[1] < b[1]
        return true
    elseif a[1] > b[1]
        return false
    else
        return a[2] < b[2]
    end
end

function central_range(min, max, n)
    dx = (max - min) / n
    range(min + 0.5dx, stop=max - 0.5dx, length=n)
end

function partition_work(N)
    k = Threads.nthreads()

    ls = range(1, stop=N, length=k+1)
    map(1:k) do i
        a = round(Int, ls[i])
        if i > 1
            a += 1
        end
        b = round(Int, ls[i+1])
        a:b
    end
end

function unique_points(points, tol)
    out = Vector{eltype(points)}()
    tol2 = tol * tol
    for p in points
        duplicate = false
        for v in out
            if sum(abs2, p .- v) < tol2
                duplicate = true
                break
            end
        end
        if !duplicate
            push!(out, p)
        end
    end
    out
end
"""
    rangeindex(range, x)

Returns an integer `i` such that `abs(range[i] - x) ≤ 0.5`
"""
rangeindex(range::StepRangeLen, x) = round(Int, (x - range[1]) / step(range)) + 1


"""
    realcoefficients(p::AbstractPolynomial)

Return the coefficient vector of `p` where each coefficient is converted to  a `Float64`.
"""
function realcoefficients(p::MP.AbstractPolynomial)::Vector{Float64}
    map(term -> float(MP.coefficient(term)), MP.terms(p))
end

"""
    mapcoefficients(f, p::AbstractPolynomial)

Construct a new polynomial from `p` where each `f` is applied on each coefficient.
"""
function mapcoefficients(f, p::MP.AbstractPolynomial)
    MP.polynomial(map(t -> f(MP.coefficient(t)) * MP.monomial(t), MP.terms(p)))
end


unzip(xs::Vector{SVector{2,T}}) where {T} = map(x -> x[1], xs), map(x -> x[2], xs)
unzip(xs::Vector{Tuple{T,T}}) where {T} = map(x -> x[1], xs), map(x -> x[2], xs)

logabs(x) = log(abs(x))


inIJulia() = isinteractive() && isdefined(Main, :IJulia) && isdefined(Main.IJulia, :clear_output)
plotsdefined() = isinteractive() && isdefined(Main, :Plots) && isdefined(Main.Plots, :plot)
inAtom() = isinteractive() && isdefined(Main, :Atom)



"""
    intersection_halfray_limits(origin::Point, d::Direction, xmin, xmax, ymin, ymax)

Compute the intersection of the halfray starting at `origin` in direction `d` with the
bounding window defined by `xmin`, `xmax`, `ymin` and `ymax`.
This returns also two (opposite) `Direction`s which are parallel to the window border
"""
function intersection_halfray_limits(origin::Point, d::Direction, xmin, xmax, ymin, ymax)
    x0, y0 = origin
    @assert xmin <= x0 <= xmax
    @assert ymin <= y0 <= ymax

    dx, dy = d

    if dx ≥ 0 && dy ≥ 0
        λx = (xmax - x0) / dx
        λy = (ymax - y0) / dy
    elseif dx < 0 && dy ≥ 0
        λx = (xmin - x0) / dx
        λy = (ymax - y0) / dy
    elseif dx ≥ 0 && dy < 0
        λx = (xmax - x0) / dx
        λy = (ymin - y0) / dy
    else
        λx = (xmin - x0) / dx
        λy = (ymin - y0) / dy
    end

    intersection = origin + min(λx, λy) * d
    # left or right border
    if λx < λy
        borderdir1 = SVector(0.0, 1.0)
        borderdir2 = -SVector(0.0, 1.0)
    # upper or lower border
    elseif λx > λy
        borderdir1 = SVector(1.0, 0.0)
        borderdir2 = -SVector(1.0, 0.0)
    # now we have to handle corners
    # upper right
    elseif dx ≥ 0 && dy ≥ 0
        borderdir1 = -SVector(1.0, 0.0)
        borderdir2 = -SVector(0.0, 1.0)
    # upper left
    elseif dx < 0 && dy ≥ 0
        borderdir1 = SVector(1.0, 0.0)
        borderdir2 = -SVector(0.0, 1.0)
    # lower right
    elseif dx ≥ 0 && dy < 0
        borderdir1 = -SVector(1.0, 0.0)
        borderdir2 = SVector(0.0, 1.0)
    # lower left
    else
        borderdir1 = SVector(1.0, 0.0)
        borderdir2 = SVector(0.0, 1.0)
    end
    intersection, borderdir1, borderdir2
end

"""
    fillrect!(bitmatrix, px, py)

Assigns `true` to the indices of `bitmatrix` which are inside the rectangle defined by `px` and `py`.
"""
function fillrect!(M::BitMatrix, px::SVector{4, Int}, py::SVector{4, Int})
    m, n = size(M)
    left, right = max(1, minimum(px)), min(n, maximum(px))
    @inbounds for x=left:right
        y1 = nothing
        y2 = nothing

        j = 4
        for i=1:4
            if (px[i] <= x <= px[j]) || (px[j] <= x <= px[i])
                # special case: adding the whole cut to ys
                if px[i] == px[j]
                    y1 = py[i]
                    y2 = py[j]
                else
                    y = round(Int, py[i] + (x - px[i]) / (px[j] - px[i]) * (py[j] - py[i]))
                    if y1 === nothing
                        y1 = y
                    elseif y2 === nothing && y != y1
                        y2 = y
                    end
                end
            end
            j = i
        end

        # if we only got a corner only add that value
        if y2 === nothing && 1 <= y1 <= m
            M[y1, x] = true
        elseif y1 !== nothing && y2 !== nothing
            a = max(min(y1, y2), 1)
            b = min(max(y1, y2), m)
            for k in a:b
                M[k, x] = true
            end
        end
    end
    return M
end

"""
    distance(p1, p2, p0)

Computes the distance of `p0` to the line through `p1` and `p2`.
"""
function distance(p1, p2, p0)
    x1, y1 = p1
    x2, y2 = p2
    x0, y0 = p0

    abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / sqrt((y2-y1)^2+(x2-x1)^2)
end


"""
    polyarea(p...)

Compute the area of the (simple) polygon defined by the vertices p.
*Note*: The last vertices does not to be the same as the first.
"""
polyarea(ps...) = polyarea(ps)
function polyarea(ps::NTuple{N, SVector{2, T}}) where {N, T}
    area = zero(T)
    for i=1:N-1
        @inbounds area += ps[i] × ps[i+1]
    end
    area += ps[N] × ps[1]
    0.5area
end

function iscontained(vertices, p)
    on_boundary(vertices, p) || in_interior(vertices, p)
end

function on_boundary(vertices, p)
    n = length(vertices)
    j = n
    i = 1
    while i ≤ n
        a, b = vertices[j], vertices[i]
        if p == a || p == b
            return true
        end
        if a[1] == b[1] == p[1] && 0 ≤ (p[2] - b[2]) / (a[2] - b[2]) ≤ 1
            return true
        elseif a[2] == b[2] == p[2] && 0 ≤ (p[1] - b[1]) / (a[1] - b[1]) ≤ 1
            return true
        end

        t1, t2 = (p - b) ./ (a - b)
        if abs(t1 - t2) < 1e-12 && 0 ≤ t1 ≤ 1
            return true
        end
        j = i
        i += 1
    end
    return false
end

function in_interior(vertices, p)
    n = length(vertices)
    j = n
    i = 1
    x, y = p
    c = false
    while i ≤ n
        vi_x, vi_y = vertices[i]
        vj_x, vj_y = vertices[j]
        if ((vi_y > y) != (vj_y > y)) && (x < (vj_x - vi_x) * (y - vi_y) / (vj_y - vi_y) + vi_x)
            c = !c
        end

        j = i
        i += 1
    end
    c
end


"""
    lineintersection(l1::NTuple{2, SVector{2, T}}, l2::NTuple{2, SVector{2, T}})

Compute the intersection between the lines `l1` and `l2`. `l1` and `l2`
are each defined by two points.
"""
function lineintersection(l1::NTuple{2, SVector{2, T}}, l2::NTuple{2, SVector{2, T}}) where T
    a1, b1 = l1
    a2, b2 = l2

    @inbounds x_a1_b1 = (a1[1] - b1[1])
    @inbounds y_a1_b1 = (a1[2] - b1[2])
    @inbounds x_a2_b2 = a2[1] - b2[1]
    @inbounds y_a2_b2 = (a2[2] - b2[2])

    a1_x_b1 = a1 × b1
    a2_x_b2 = a2 × b2

    denom = inv(x_a1_b1 * y_a2_b2 - y_a1_b1 * x_a2_b2)

    x = (a1_x_b1 * x_a2_b2 - x_a1_b1 * a2_x_b2) * denom
    y = (a1_x_b1 * y_a2_b2 - y_a1_b1 * a2_x_b2) * denom

    SVector{2}(x, y)
end


function equal_aspect_ratio(xmin, xmax, ymin, ymax)
    Δy = ymax - ymin
    Δx = xmax - xmin

    ratio = Δy / Δx

    if ratio ≈ 1.0
        return (xmin, xmax, ymin, ymax)
    elseif ratio > 1.0
        # y is longer than x -> increase x
        δ = 0.5 * (Δy - Δx)
        xmin -= δ
        xmax += δ

        return (xmin, xmax, ymin, ymax)
    else
        δ = 0.5 * (Δy - Δx)
        ymin -= δ
        ymax += δ

        return (xmin, xmax, ymin, ymax)
    end
end
