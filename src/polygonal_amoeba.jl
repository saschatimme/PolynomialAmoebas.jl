export domain, accuracy, spine

"""
    domain(A::PolygonalAmoeba)

The domain `(xmin, xmax, ymin, ymax)` for which the amoeba is defined.
"""
domain(A::PolygonalAmoeba) = A.domain

"""
    accuracy(A::PolygonalAmoeba)

The accuracy  of the approximation.
"""
accuracy(A::PolygonalAmoeba) = A.error_bound

"""
    spine(A::PolygonalAmoeba)

The spine which was used to compute `A`.
"""
spine(A::PolygonalAmoeba) = A.spine


function Base.show(io::IO, A::PolygonalAmoeba)
    println(io, "PolygonalAmoeba:")
    println(io, " * accuracy: ", accuracy(A))
end

function Base.show(io::IO, mime::MIME"text/html", A::PolygonalAmoeba)
    if inIJulia() && plotsdefined()
        show(io, mime, Main.Plots.plot(A))
    else
        show(io, A)
    end
end


@recipe function plot(A::PolygonalAmoeba; spine=false)
    xmin, xmax, ymin, ymax = A.domain
    xlims --> (xmin, xmax)
    ylims --> (ymin, ymax)
    aspect_ratio --> :equal
    for (p1, p2, p3, p4) in A.polygons
        @series begin
            fillalpha--> get(plotattributes, :seriesalpha, 1.0)
            linecolor --> :dodgerblue
            fillcolor --> :dodgerblue
            linewidth --> 1.0
            linealpha --> get(plotattributes, :seriesalpha, get(plotattributes, :fillalpha, 1.0))
            legend --> false
            st := :shape
            x1, y1 = p1
            x2, y2 = p2
            x3, y3 = p3
            x4, y4 = p4
            x5, y5 = p1

            [x1, x2, x3, x4, x5], [y1, y2, y3, y4, y5]
        end
    end
    if spine
        @series begin
            A.spine
        end
    end
end
