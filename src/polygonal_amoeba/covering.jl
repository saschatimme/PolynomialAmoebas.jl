"""
    CoverPoly(p1::Pillar, p2::Pillar)

A polygon defined by the two pillars `p1` and `p2` where `p1` and `p2`
point to points ``x ∈ E_α(f)``.
"""
struct CoverPoly
    p1::Pillar
    p2::Pillar
end

struct CoveringPolygon
    p1_id::Int
    p2_id::Int
    order::SVector{2, Int}
    left_id::Int
    right_id::Int
end

mutable struct Covering
    pillars::Vector{Pillar}
    polygons::Dict{Int, CoveringPolygon}
    id_counter::Int # id of latest inserted polygon, increment this to get unique ids
end
function Covering(pillars::Vector{Pillar}, polygons::Vector{CoveringPolygon})
    id_counter = length(polygons)
    Covering(pillars, Dict(enumerate(polygons)), id_counter)
end

Covering(S::Spine2D; kwargs...) = Covering(S, S.amoeba_approximation.grid; kwargs...)
Covering(S::Spine2D, grid::Grid2D; kwargs...) = Covering(S, limits(grid)...; kwargs...)
function Covering(S::Spine2D, xmin::Real, xmax::Real, ymin::Real, ymax::Real)
    handle_halfray_intersections=true
    spine = S.curve

    # we have to ensue that the points of the complement complements are inside the
    # window. So we enlarge it if necessary
    CCs = components_complement(S)
    lims = covering_limits(CCs, xmin, xmax, ymin, ymax)
    intersect_vertices, vertices, linesegments = vertices_linesegments(S.curve, lims)

    order_pillars_dict = make_order_pillars_dict(CCs, vertices, linesegments)
    sorted_order_pillars = sort_order_pillars(order_pillars_dict, CCs, intersect_vertices)

    # We need to take care of the extra covering at the intersections
    # due to our sorting these are origins of the first and the last pillar
    if handle_halfray_intersections
        for (cc, pillars) in zip(CCs, sorted_order_pillars)
            if cc.bounded
                continue
            end
            p_start = intersection_cover(pillars[1], pillars[2], xmin, xmax, ymin, ymax)
            pushfirst!(pillars, p_start)
            p_end = intersection_cover(pillars[end], pillars[end-1], xmin, xmax, ymin, ymax)
            push!(pillars, p_end)
        end
    end
    # sorted_order_pillars is now a vector where each entry is a vector containing
    # all pillars to a given order (sorted in positive direction).
    # We can now construct the polygonal cover

    k = 0
    pillar_vec = Vector{Pillar}()
    polygons = Vector{CoveringPolygon}()
    poly_id = 0
    pillar_id = 0
    for (cc, pillars) in zip(CCs, sorted_order_pillars)
        n = length(pillars)
        if cc.bounded
            j = n
            i = 1
            while i ≤ n
                p1, p2 = pillar_id + j, pillar_id + i
                left_id = poly_id + j
                if i == n
                    right_id = poly_id + 1
                else
                    right_id = poly_id + i + 1
                end
                poly = CoveringPolygon(p1, p2, cc.order, left_id, right_id)
                push!(polygons, poly)
                j = i
                i += 1
            end
            poly_id += n
        else
            for i=1:n-1
                p1, p2 = pillar_id + i, pillar_id + i + 1
                left_id = i == 1 ? 0 : poly_id + i - 1
                right_id = i == n - 1 ? 0 : poly_id + i + 1
                poly = CoveringPolygon(p1, p2, cc.order, left_id, right_id)
                push!(polygons, poly)
            end
            poly_id += n - 1
        end
        pillar_id += n
        append!(pillar_vec, pillars)
    end
    Covering(pillar_vec, polygons)
end


function covering_limits(CCs, xmin, xmax, ymin, ymax)
    x_complement_points = first.(point.(CCs))
    y_complement_points = last.(point.(CCs))
    xmin = min(xmin, minimum(x_complement_points))
    xmax = max(xmax, maximum(x_complement_points))
    ymin = min(ymin, minimum(y_complement_points))
    ymax = max(ymax, maximum(y_complement_points))

    return xmin, xmax, ymin, ymax
end

function vertices_linesegments(spine, limits)
    linesegments = Vector{SVector{2, Point}}()
    all_vertices = Vector{SVector{2,Float64}}()
    xmin, xmax, ymin, ymax = limits
    for (origin, direction) in halfrays(spine)
        intersec, _, _ = intersection_halfray_limits(origin, direction, xmin, xmax, ymin, ymax)
        push!(linesegments, SVector(origin, intersec))
        push!(all_vertices, intersec)
    end

    intersect_vertices = copy(all_vertices)
    append!(linesegments, segments(spine))
    append!(all_vertices, vertices(spine))

    intersect_vertices, all_vertices, linesegments
end
function make_order_pillars_dict(CCs, vertices, linesegments)
    order_pillars_dict = map(cc -> (cc.order, Vector{Pillar}()), CCs) |> Dict
    for a in vertices
        for cc in CCs
            p = cc.x
            Δap = p - a
            nointersection = all(linesegments) do l
                !linesegments_intersect_pp_pd(l, (a, Δap))
            end
            # we found a proper segment
            if nointersection
                push!(order_pillars_dict[cc.order], (a, Δap))
            end
        end
    end
    return order_pillars_dict
end

function update_neighbor(p::CoveringPolygon, nb_id::Int, left_or_right::Symbol)
    if left_or_right == :left
        CoveringPolygon(p.p1_id, p.p2_id, p.order, nb_id, p.right_id)
    elseif left_or_right == :right
        CoveringPolygon(p.p1_id, p.p2_id, p.order, p.left_id, nb_id)
    else
        throw(error("Expected `:left` or `:right`, got `:$(left_or_right)`"))
    end
end

function cartesian_polygon(covering::Covering, cp::CoveringPolygon)
    cartesian_polygon(cp, covering.pillars)
end
function cartesian_polygon(cp::CoveringPolygon, pillars::Vector{Pillar})
    p1, p2 = pillars[cp.p1_id], pillars[cp.p2_id]
    (p2[1], sum(p2), sum(p1), p1[1])
end

area(covering::Covering, cp::CoveringPolygon) = area(cp, covering.pillars)
area(cp::CoveringPolygon, pillars) = polyarea(cartesian_polygon(cp, pillars))

polygons(covering::Covering) = covering.polygons
pillar(covering::Covering, id::Int) = covering.pillars[id]
polygon(covering::Covering, id::Int) = covering.polygons[id]

function _set!(covering, id::Int, pillar::Pillar)
    covering.pillars[id] = pillar
    covering
end
function _set!(covering, id::Int, poly::CoveringPolygon)
    covering.polygons[id] = poly
    covering
end

"""
    new_poly_id!(covering)

Create a new poly id.
"""
function new_poly_id!(covering::Covering)
    covering.id_counter += 1
    covering.id_counter
end
"""
    add_pillar!(covering, pillar)

Add a new pillar to the covering and return its id.
"""
function add_pillar!(covering::Covering, pillar::Pillar)
    push!(covering.pillars, pillar)
    length(covering.pillars)
end

"""
    add_polygon!(covering, id, poly)

Add a polygon with the given `id` to the covering.
"""
function add_polygon!(covering::Covering, id::Int, poly::CoveringPolygon)
    _set!(covering, id, poly)
    covering
end

"""
    remove_polygon!(covering, id, poly)

Remove a polygon with the given `id` from the covering.
"""
function remove_polygon!(covering::Covering, id::Int)
    delete!(covering.polygons, id)
    covering
end

"""
    update_polygon_neighbor!(covering, id, new_neighbor_id, left_right)

Update the left or right neighbor of the polygon with the given `id` to be
`new_neighbor_id`.
"""
function update_polygon_neighbor!(covering::Covering, id::Int, new_neighbor_id::Int, left_right::Symbol)
    poly = covering.polygons[id]
    covering.polygons[id] = update_neighbor(poly, new_neighbor_id, left_right)
    covering
end

function cartesian_polygons(covering::Covering)
    map(cp -> cartesian_polygon(cp, covering.pillars), values(covering.polygons))
end


function sort_order_pillars(order_pillars_dict, CCs, intersect_vertices)
    map(CCs) do cc
        pillars = order_pillars_dict[cc.order]
        vertices = first.(pillars)

        hull_indices = convexhull_vertices(vertices)
        hull = vertices[hull_indices]
        dirs = (last.(pillars))[hull_indices]
        sorted_pillars = collect(zip(hull, dirs))
        if cc.bounded
            return sorted_pillars
        else
            n = length(hull)
            if any(v -> v == hull[n], intersect_vertices) &&
               any(v -> v == hull[1], intersect_vertices)
               return sorted_pillars
            end

            for i=1:n-1
                if any(v -> v == hull[i], intersect_vertices) &&
                   any(v -> v == hull[i+1], intersect_vertices)
                   return [sorted_pillars[i+1:end]; sorted_pillars[1:i]]
                end
            end
        end
    end
end

function intersection_cover(intersection_pillar, node_pillar, xmin, xmax, ymin, ymax)
    intersection, Δintersection_p = intersection_pillar
    node, Δnode_p = node_pillar
    direction = intersection - node
    # we recompute the intersection to get the canidates for the pillar direciton
    _, d1, d2 = intersection_halfray_limits(node, 2direction, xmin, xmax, ymin, ymax)

    # we need to have a proper line segment to apply `pillars_onsameside`
    if pillars_onsameside((node, Δnode_p), (intersection, d1))
        dir = d1
    else
        dir = d2
    end
    # d is the projection of `dir` onto the line parallel to `direction`
    # through `p`
    d = (Δintersection_p × direction) / (dir × direction) * dir
    # # we maybe have to scale d, this can happen for
    # λ = scale_to_fit(intersection, d, xmin, xmax, ymin, ymax)
    (intersection, d)
end

"""
    scalepillars!(f, covering::Covering)

Scale each pillar of `covering` by a factor `t` returned by the function `f`.
`f` has as an input the origin and direction of the pillar.
"""
function scalepillars!(f, covering::Covering)
    N = length(covering.pillars)
    for k=1:N
        origin, dir = pillar(covering, k)
        t = f(origin, dir)
        covering.pillars[k] = (origin, t * dir)
    end
    covering
end



function split!(covering::Covering, polygon_id::Int, new_pillar::Pillar)
    pillar_id = add_pillar!(covering, new_pillar)
    poly = polygon(covering, polygon_id)

    new_poly_id_left = new_poly_id!(covering)
    new_poly_id_right = new_poly_id!(covering)

    # create new polygons
    left = CoveringPolygon(
        poly.p1_id, pillar_id, poly.order, poly.left_id, new_poly_id_right)
    right = CoveringPolygon(
        pillar_id, poly.p2_id, poly.order, new_poly_id_left, poly.right_id)

    # now we have to update the left and right neighbours of p
    if poly.left_id != 0
        update_polygon_neighbor!(covering, poly.left_id, new_poly_id_left, :right)
    end
    if poly.right_id != 0
        update_polygon_neighbor!(covering, poly.right_id, new_poly_id_right, :left)
    end

    remove_polygon!(covering, polygon_id)
    add_polygon!(covering, new_poly_id_left, left)
    add_polygon!(covering, new_poly_id_right, right)

    (new_poly_id_left, new_poly_id_right)
end


# AREA
function intersection_point(covering::Covering, poly::CoveringPolygon, ε)
    l1, l2 = approximation_lines(covering, poly, ε)
    lineintersection(l1, l2)
end
function intersection_point(c::Covering, poly_id::Int, ε)
    intersection_point(c, polygon(covering, poly_id), ε)
end

function approximation_lines(covering::Covering, poly::CoveringPolygon, ε)
    if poly.left_id != 0 && poly.right_id != 0
        left = polygon(covering, poly.left_id)
        right = polygon(covering, poly.right_id)

        # we have to account for the possible ε error
        a1, v1 = pillar(covering, poly.p1_id)
        a2, v2 = pillar(covering, poly.p2_id)

        b1 = a1 + (1 - ε / norm(v1)) * v1
        b2 = a2 + (1 - ε / norm(v2)) * v2

        l1 = (b1, sum(pillar(covering, left.p1_id)))
        l2 = (b2, sum(pillar(covering, right.p2_id)))

        return l1, l2
    elseif poly.left_id == 0
        # We have one of these special polygons
        right = polygon(covering, poly.right_id)

        # a1 == a2 and (a1 + t*v1) is part of the boundary
        a1, v1 = pillar(covering, poly.p1_id)
        a2, v2 = pillar(covering, poly.p2_id)

        b2 = a2 + (1 - ε / norm(v2)) * v2

        l1 = (a1, a1 + v1)
        l2 = (b2, sum(pillar(covering, right.p2_id)))

        return l1, l2
    else
        #poly.right_id == 0
        # We have one of these special polygons
        left = polygon(covering, poly.left_id)

        # a1 == a2 and (a2 + t*v2) is part of the boundary
        a1, v1 = pillar(covering, poly.p1_id)
        a2, v2 = pillar(covering, poly.p2_id)

        b1 = a1 + (1 - ε / norm(v1)) * v1

        l2 = (a2, a2 + v2)
        l1 = (b1, sum(pillar(covering, left.p1_id)))

        return l1, l2
    end
end



function approximated_local_error(covering::Covering, poly_id::Int, ε)
    approximated_local_error(covering, polygon(covering, poly_id), ε)
end
function approximated_local_error(covering::Covering, p::CoveringPolygon, ε)
    z = intersection_point(covering, p, ε)
    poly = cartesian_polygon(covering, p)
    if iscontained(poly, z)
        if p.left_id == 0
            x1 = sum(pillar(covering, p.p1_id))
            a2, v2 = pillar(covering, p.p2_id)
            y2 = a2 + (1 - ε / norm(v2)) * v2
            x2 = a2 + v2
            polyarea(z, y2, x2, x1)
        elseif p.right_id == 0
            a1, v1 = pillar(covering, p.p1_id)
            x1 = a1 + v1
            x2 = sum(pillar(covering, p.p2_id))

            y1 = a1 + (1 - ε / norm(v1)) * v1

            polyarea(z, x2, x1, y1)
        else
            a1, v1 = pillar(covering, p.p1_id)
            a2, v2 = pillar(covering, p.p2_id)
            y1 = a1 + (1 - ε / norm(v1)) * v1
            y2 = a2 + (1 - ε / norm(v2)) * v2
            x1, x2 = a1 + v1, a2 + v2
            return polyarea(z, y2, x2, x1, y1)
        end
    else
        return polyarea(poly)
    end
end

"""
    side(A::Point, B::Point, p::Point)

`side(A,B,p)>0` if `p` lies on one side of the line through `A` and `B` and `side(A,B,p)<0`
if `p` lies on te other side and `side(A,B,p) == 0` if `p` lies on the line.
"""
side(A, B, p) = (p[1] - A[1]) * (B[2] - A[2]) - (p[2] - A[2]) * (B[1] - A[1])

"""
    pillars_onsameside(p1::Pillar, p2::Pillar)

Check if two pillars point to the same side as defined by the line between the origins of
`p1` and `p2`.
"""
function pillars_onsameside(p1::Pillar, p2::Pillar)
     side(p1[1], p2[1], p1[1] + p1[2]) * side(p1[1], p2[1], p2[1] + p2[2]) >= 0
end

"""
    scale_to_fit(origin, d, xmin, xmax, ymin, ymax)

Return a scaling factor ``λ`` such that `origin + λ * d` is on the boundary of the window
determined by `xmin`, `xmax`, `ymin` and `ymax`.
"""
function scale_to_fit(origin::Point, d::Direction, xmin, xmax, ymin, ymax)
        x0, y0 = origin
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

        min(λx, λy)
end


"""
    linesegments_intersect(q::Point, s::Direction, p::Point, r::Direction)

Determines whether the line segments `(q, q + s)` and `(p, p + r)` intersect
"""
function linesegments_intersect(q::Point, s::Direction, p::Point, r::Direction)
    # algorithm taken from
    # https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
    unum = (q - p) × r
    tnum = (q - p) × s
    denom = r × s

    # colinear
    if unum ≈ 0.0 && denom ≈ 0
        norm2 = r ⋅ r
        t0 = (q - p) ⋅ r / norm2
        t1 = t0 + s ⋅ r / norm2
        if s ⋅ r < 0
            return 0 ≤ t0 ≤ 1 || 0 ≤ t1 ≤ 1 || t0 ≤ 0 && 1 ≤ t1
        else
            return 0 ≤ t0 ≤ 1 || 0 ≤ t1 ≤ 1 || t1 ≤ 0 && 1 ≤ t0
        end
    elseif denom ≈ 0
        return false
    else
        t = tnum / denom
        u = unum / denom

        # we only count proper intersections, not if endpoints coincide
        return 0 < t < 1 && 0 < u < 1
    end
end
@inline function linesegments_intersect(a::Tuple{Point, Direction}, b::Tuple{Point, Direction})
    linesegments_intersect(a[1], a[2], b[1], b[2])
end

function linesegments_intersect_pp_pd(a::SVector{2, Point}, b::Tuple{Point, Direction})
    linesegments_intersect(a[1], a[2] - a[1], b[1], b[2])
end


@recipe function plot(covering::Covering)
    fillalpha--> 0.4
    linecolor --> :dodgerblue
    fillcolor --> :dodgerblue
    linewidth --> 1.0
    legend --> false
    for (p1, p2, p3, p4) in cartesian_polygons(covering)
        @series begin
            st := :shape
            x1, y1 = p1
            x2, y2 = p2
            x3, y3 = p3
            x4, y4 = p4
            x5, y5 = p1

            [x1, x2, x3, x4, x5], [y1, y2, y3, y4, y5]
        end
    end
end
