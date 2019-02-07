export Grid3D

"""
    Grid3D(;xlims=(-1,1), ylim=(-1, 1), zlims=(-1, 1), res=(21, 21, 21))

Create a 3D grid.
"""
function Grid3D(;xlims=(-1, 1), ylims=(-1, 1), zlims=(-1, 1), res=(21, 21, 21))
    Grid3D(
        central_range(xlims..., res[1]),
        central_range(ylims..., res[2]),
        central_range(zlims..., res[3]))
end

function Grid3D(lims, res::Int)
    Grid3D(
        central_range(lims..., res),
        central_range(lims..., res),
        central_range(lims..., res))
end

Base.broadcastable(G::Grid3D) = Ref(G)
Base.size(G::Grid3D) = (length(G.xrange), length(G.yrange), length(G.zrange))
Base.size(G::Grid3D, i) = size(G)[i]
Base.length(G::Grid3D) = prod(size(G))
Base.ndims(G::Grid3D) = 3

@inline function Base.getindex(G::Grid3D, s::Integer)
    i, j, k = _ind2sub(size(G), s)
    G[i, j, k]
end
Base.getindex(G::Grid3D, I::CartesianIndex{3}) = G[Tuple(I)...]

@inline function Base.getindex(G::Grid3D, i::Integer, j::Integer, k::Integer)
    G.xrange[i], G.yrange[j], G.zrange[k]
end
Base.getindex(G::Grid3D, ijk::NTuple{3, Int}) = G[ijk[1], ijk[2], ijk[3]]
Base.getindex(G::Grid3D, xyz::NTuple{3, Float64}) = gridpoint(G, xyz[1], xyz[2], xyz[3])

xcoordinates(G::Grid3D) = G.xrange
ycoordinates(G::Grid3D) = G.yrange
zcoordinates(G::Grid3D) = G.zrange

xcoordinate(G::Grid3D, i) = G.xrange[i]
ycoordinate(G::Grid3D, i) = G.yrange[i]
zcoordinate(G::Grid3D, i) = G.urange[i]

xgridpoint(G::Grid3D, x::Float64) = rangeindex(G.xrange, x)
ygridpoint(G::Grid3D, y::Float64) = rangeindex(G.yrange, y)
zgridpoint(G::Grid3D, z::Float64) = rangeindex(G.zrange, z)
gridpoint(G::Grid3D, x, y, z) = (xgridpoint(G, x), ygridpoint(G, y), zgridpoint(G, z))
