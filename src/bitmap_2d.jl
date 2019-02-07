"""
    Bitmap2D(grid)

Create an empty `Bitmap2D` from an grid.
"""
Bitmap2D(grid::Grid2D) = Bitmap2D(falses(size(grid)...), grid)

Base.copy(M::Bitmap2D) = Bitmap2D(copy(M.data), M.grid)
Base.length(M::Bitmap2D) = length(M.grid)
Base.size(M::Bitmap2D) = size(M.grid)
Base.size(M::Bitmap2D, i::Int64) = Base.size(M.grid, i)
Base.ndims(M::Bitmap2D) = 2
Base.getindex(M::Bitmap2D, i) = getindex(M.data, i)
Base.getindex(M::Bitmap2D, i, j) = getindex(M.data, i, j)
Base.setindex!(M::Bitmap2D, x::Bool, i) = setindex!(M.data, x, i)
Base.setindex!(M::Bitmap2D, x::Bool, i, j) = setindex!(M.data, x, i, j)

Base.intersect(M₁::Bitmap2D, M₂::Bitmap2D) = Bitmap2D(M₁.data .& M₂.data, M₂.grid)

function safe_getindex(M::Bitmap2D, i, j)
    if 1 ≤ i ≤ size(M.data, 1) && 1 ≤ j ≤ size(M.data, 2)
        return getindex(M.data, i, j)
    end
    return false
end

function safe_setindex!(M::Bitmap2D, x::Bool, i, j)
    if 1 ≤ i ≤ size(M.data, 1) && 1 ≤ j ≤ size(M.data, 2)
        setindex!(M, x, i, j)
    end
end

@recipe function plot(M::Bitmap2D; transparent=true, scatter=false)
    seriescolor --> :dodgerblue

    if scatter
        xs, ys = _scatterdata(M)

        @series begin
            st = :scatter
            xs, ys
        end

    else
        clr = Colors.parse(Colors.Colorant, plotattributes[:seriescolor])
        cmap = PlotUtils.ColorGradient([clr, clr])
        fillcolor := cmap
        colorbar := false
        aspect_ratio --> :equal
        data = map(x -> x ? 1.0 : NaN, M.data)
        if !transparent
            cmap = PlotUtils.ColorGradient([Colors.parse(Colors.Colorant, :white), clr])
            fillcolor := cmap
            data = map(x -> x ? 1.0 : 0.0, M.data)
        end

        (xmin, xmax, ymin, ymax) = limits(M.grid)

        @series begin
            st := :heatmap
            grid  --> false
            xlims --> (xmin, xmax)
            ylims --> (ymin, ymax)
            xcoordinates(M.grid), ycoordinates(M.grid), data'
        end
    end
end

function _scatterdata(D::Bitmap2D)
    xs = Float64[]
    ys = Float64[]

    for l in eachindex(D.data)
        if D[l]
            i, j = ind2sub(D.data, l)
            push!(xs, D.grid.xrange[i])
            push!(ys, D.grid.yrange[j])
        end
    end
    (xs, ys)
end

fillrect!(M::Bitmap2D, px, py) = fillrect!(M.data, px, py)


#
# function Base.show(io::IO, M::Bitmap2D)
#     #println(typeof(A), ":")
#     println(io, "$(size(M, 1))×$(size(M, 2)) Bitmap2D with window [$(M.xmin), $(M.xmax), $(M.ymin), $(M.ymax)]")
# end
