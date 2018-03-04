"""
    Bitmap3D(grid::Grid3d)

Create an empty `Bitmap3D` from an grid.
"""
Bitmap3D(grid::Grid3D) = Bitmap3D(falses(size(grid)...), grid)

Base.copy(M::Bitmap3D) = Bitmap3D(copy(M.data), grid)

Base.size(M::Bitmap3D) = size(M.grid)
Base.size(M::Bitmap3D, i::Int64) = Base.size(M.grid, i)
Base.ndims(M::Bitmap3D) = 3
Base.getindex(M::Bitmap3D, i) = getindex(M.data, i)
Base.getindex(M::Bitmap3D, i, j, k) = getindex(M.data, i, j, k)
Base.setindex!(M::Bitmap3D, x::Bool, i) = setindex!(M.data, x, i)
Base.setindex!(M::Bitmap3D, x::Bool, i, j, k) = setindex!(M.data, x, i, j, k)
Base.eachindex(M::Bitmap3D) = eachindex(M.data)

function safe_getindex(M::Bitmap3D, i, j, k)
    if 1 ≤ i ≤ size(M.data, 1) &&
       1 ≤ j ≤ size(M.data, 2) &&
       1 ≤ k ≤ size(M.data, 3)
        return getindex(M, i, j, k)
    end
    return false
end
function safe_setindex!(M::Bitmap3D, x::Bool, i, j, k)
    if 1 ≤ i ≤ size(M.data, 1) &&
       1 ≤ j ≤ size(M.data, 2) &&
       1 ≤ k ≤ size(M.data, 3)
        setindex!(M, x, i, j, k)
    end
end


@recipe function plot(M::Bitmap3D; boundarycolor=:dodgerblue)
    inner, boundary = _scatterdata(M)

    @series begin
        seriescolor --> :dodgerblue
        legend --> false
        st := :scatter3d
        inner
    end

    @series begin
        seriescolor --> boundarycolor
        legend --> false
        st := :scatter3d
        boundary
    end
end

function _scatterdata(D::Bitmap3D)
    xs = Float64[]
    ys = Float64[]
    zs = Float64[]

    bxs = Float64[]
    bys = Float64[]
    bzs = Float64[]
    for l in eachindex(D.data)
        if D[l]
            if !isboundary_point(D.data, l)
                i, j, k = ind2sub(D.data, l)
                push!(xs, D.grid.xrange[i])
                push!(ys, D.grid.yrange[j])
                push!(zs, D.grid.zrange[k])
            else
                i, j, k = ind2sub(D.data, l)
                push!(bxs, D.grid.xrange[i])
                push!(bys, D.grid.yrange[j])
                push!(bzs, D.grid.zrange[k])
            end
        end
    end
    (xs, ys, zs), (bxs, bys, bzs)
end

function isboundary_point(data, l)
    i, j, k = ind2sub(data, l)
    i0 = max(i - 1, 1)
    i1 = min(i + 1, size(data, 1))

    j0 = max(j - 1, 1)
    j1 = min(j + 1, size(data, 2))

    k0 = max(k - 1, 1)
    k1 = min(k + 1, size(data, 3))
    for ii = i0:i1, jj = j0:j1, kk=k0:k1
        if !data[ii, jj, kk]
            return true
        end
    end
    return false
end
