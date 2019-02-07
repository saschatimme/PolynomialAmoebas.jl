export ConvexHull, convexhull

const spatial = PyNULL()

function __init__()
    copy!(spatial, pyimport_conda("scipy.spatial", "scipy"))
end

struct ConvexHull{T<:Real}
    points::Matrix{T}
    vertices::Vector{Int}
    simplices::Vector{Vector{Int}}
    equations::Matrix{T}
end

addone(xs) = map(x -> x + 1, xs)

function convexhull(x::Matrix{T}) where T<:Real
    py = spatial[:ConvexHull](x)
    points = convert(Matrix{T}, py["points"])
    vertices = addone(convert(Vector{Int}, py["vertices"]))
    simplices = addone.(convert(Vector{Vector{Int}}, py["simplices"]))
    equations = convert(Matrix{T}, py["equations"])
    ConvexHull(points, vertices, simplices, equations)
end
"""
    convexhull_vertices(line::Vector{Tuple{Int, Int}})

Returns a list of indices of the vertices of the convex hull.
"""
function convexhull_vertices(line::Vector{<:Union{SVector{2, T}, NTuple{2, T}}}) where T<:Real
    m = length(line)
    M = Matrix{T}(undef, m, 2)
    for k = 1:m
        i, j = line[k]
        M[k, 1] = i
        M[k, 2] = j
    end
    py = spatial[:ConvexHull](M)
    vertices = addone(convert(Vector{Int}, py["vertices"]))
    return vertices
end
