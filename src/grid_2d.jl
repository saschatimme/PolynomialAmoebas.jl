export Grid2D

"""
    Grid2D(;xlims=(-1,1), ylim=(-1, 1), res=(21, 21))

Create a 2D grid.
"""
function Grid2D(;xlims=(-1, 1), ylims=(-1, 1), res::NTuple{2, Int}=(21, 21))
    Grid2D(
        central_range(xlims..., res[1]),
        central_range(ylims..., res[2]))
end
function Grid2D(xmin, xmax, ymin, ymax, xres, yres)
    Grid2D(
        central_range(xmin, xmax, xres),
        central_range(ymin, ymax, yres))
end
# Grid2D(;lims=(-1, 1), res::Int=21) = Grid2D(lims, res)
Grid2D(lims::Tuple{<:Real, <:Real}, res::Int) = Grid2D( central_range(lims..., res), central_range(lims..., res))
function Grid2D(f::MP.AbstractPolynomial, res=600; factor=1.5)
    Grid2D(amoeba_carcase_domain_heuristic(f, factor=factor)..., res, res)
end

Base.broadcastable(G::Grid2D) = Ref(G)
Base.size(G::Grid2D) = (length(G.xrange), length(G.yrange))
Base.size(G::Grid2D, i) = size(G)[i]
Base.length(G::Grid2D) = prod(size(G))
Base.ndims(G::Grid2D) = 3
@inline function Base.getindex(G::Grid2D, k::Integer)
    i, j = _ind2sub(size(G), k)
    G[i, j]
end
Base.getindex(G::Grid2D, I::CartesianIndex{2}) = G[Tuple(I)...]
@inline function Base.getindex(G::Grid2D, i::Integer, j::Integer)
    G.xrange[i], G.yrange[j]
end
Base.getindex(G::Grid2D, ij::Tuple{<:Integer, <:Integer}) = G[ij[1], ij[2]]
Base.getindex(G::Grid2D, ij::SVector{2, <:Integer}) = G[ij[1], ij[2]]
Base.getindex(G::Grid2D, xy::NTuple{2, Float64}) = gridpoint(G, xy[1], xy[2])

xcoordinates(G::Grid2D) = G.xrange
ycoordinates(G::Grid2D) = G.yrange
xcoordinate(G::Grid2D, i) = G.xrange[i]
ycoordinate(G::Grid2D, i) = G.yrange[i]

xgridpoint(G::Grid2D, x::Float64) = rangeindex(G.xrange, x)
ygridpoint(G::Grid2D, y::Float64) = rangeindex(G.yrange, y)
gridpoint(G::Grid2D, x, y) = (xgridpoint(G, x), ygridpoint(G, y))
function limits(G::Grid2D)
    xcoords, ycoords = xcoordinates(G), ycoordinates(G)
    xmin = xcoords[1] - 0.5 * step(xcoords)
    xmax = xcoords[end] + 0.5 * step(xcoords)
    ymin = ycoords[1] - 0.5 * step(ycoords)
    ymax = ycoords[end] + 0.5 * step(ycoords)

    xmin, xmax, ymin, ymax
end

"""
    line(f, G, p0, p1)

Map the line between `p0` and `p1` on the grid `G` and call function `f` for each
grid point on that line.
This is an implemenation of
[Bresenhams's line algorithm](https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm).
"""
function line(f::Function, G::Grid2D, p0, p1)
    x, y = p0
    x1, y1 = p1
    m, n = size(G)

    dx = round(Int, abs(x1-x))
    dy = round(Int, abs(y1-y))

    sx = x < x1 ? 1 : -1
    sy = y < y1 ? 1 : -1;

    err = 0.5 * (dx > dy ? dx : -dy)

    while x != x1 || y != y1
        if 1 ≤ x ≤ m && 1 ≤ y ≤ n
            f(x, y)
        end

        e2 = err;
        if e2 > -dx
            err -= dy
            x += sx
        end
        if e2 < dy
            err += dx
            y += sy
        end
    end
end
