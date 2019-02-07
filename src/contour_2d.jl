export Contour2D, points

"""
    Contour2D

Represents a contour of a two-dimensional amoeba.
"""
struct Contour2D
    ps::Vector{SVector{2, Float64}}
    grid::Grid2D
end

"""
    points(C::Contour2D)

The sample points describing `C`.
"""
points(C::Contour2D) = C.ps

function Base.show(io::IO, A::Contour2D)
    println(io, "Contour2D:")
    println(io, " * domain: ", limits(A.grid))
    println(io, " * number of points: ", length(points(A)))
end

function Base.show(io::IO, mime::MIME"text/html", A::Contour2D)
    if inIJulia() && plotsdefined()
        show(io, mime, Main.Plots.plot(A))
    else
        show(io, A)
    end
end



@recipe function plot(C::Contour2D)
    xmin, xmax, ymin, ymax = limits(C.grid)
    @series begin

        if !haskey(plotattributes, :xlims)
            xlims --> (xmin, xmax)
        end
        if !haskey(plotattributes, :ylims)
            ylims --> (ymin, ymax)
        end
        aspect_ratio --> :equal
        # grid --> false
        legend := false
        st := :scatter
        #markershape := :rect
        markercolor --> Colors.colorant"#884EA0"
        markerstrokewidth := 0.0
        markersize --> 2.0

        xs, ys = unzip(C.ps)
        xs, ys
    end

end


function _contour2D(F::ContourFiber2D, grid::Grid2D,
    samples_off_axis,
    options::MembershipTestOptions)

    ntries = options.ntries
    iterations = options.maxiters
    tol = options.tol

    fix_axis!(F, :w1)
    Fs = [deepcopy(F) for _=2:Threads.nthreads()]
    push!(Fs, F)

    threads_ps = [Vector{SVector{2, Float64}}() for _=1:Threads.nthreads()]
    threads_ws = [Vector{Float64}() for _=1:Threads.nthreads()]

    xmin, xmax, ymin, ymax = limits(grid)
    # first pass
    w1s = xcoordinates(grid)
    w2s = range(ymin, stop=ymax, length=samples_off_axis)

    force_exit_y(x) = x[3] > ymax || x[3] < ymin

    Threads.@threads for k=1:length(w1s)
        G = Fs[Threads.threadid()]
        w1_ys = threads_ws[Threads.threadid()]
        ps = threads_ps[Threads.threadid()]
        w1 = w1s[k]
        update_fiber!(G, w1)

        for w2 in w2s
            k = 0
            while k < ntries
                k += 1
                x = SVector(2π*rand(), 2π*rand(), w2)
                res = newton(G, x, iterations, tol, force_exit_y)
                if res.converged
                    y = res.solution[3]
                    p = SVector(w1, y)
                    if isunique(w1_ys, y, step(w1s) * 0.5)
                        push!(ps, p)
                        push!(w1_ys, y)
                    end
                end
            end
        end
        empty!(w1_ys)
    end

    # second pass
    for F in Fs
        fix_axis!(F, :w2)
    end
    w1s = range(xcoordinates(grid)[1], stop=xcoordinates(grid)[end], length=samples_off_axis)
    w2s = ycoordinates(grid)
    # w2_xs = Vector{Float64}()
    force_exit_x(x) = x[3] > xmax || x[3] < xmin

    Threads.@threads for k=1:length(w2s)
        G = Fs[Threads.threadid()]
        w2_xs = threads_ws[Threads.threadid()]
        ps = threads_ps[Threads.threadid()]
        w2 = w2s[k]
        update_fiber!(G, w2)
        for w1 in w1s
            k = 0
            while k < ntries
                k += 1
                x = SVector(2π*rand(), 2π*rand(), w1)
                res = newton(G, x, iterations, tol, force_exit_x)
                if res.converged
                    x = res.solution[3]
                    p = SVector(x, w2)
                    if isunique(w2_xs, x, step(w2s) * 0.5)
                        push!(ps, p)
                        push!(w2_xs, x)
                    end
                end
            end
        end
        empty!(w2_xs)
    end

    ps = Vector{SVector{2, Float64}}()
    for i=1:Threads.nthreads()
        append!(ps, threads_ps[i])
    end
    tol = 0.5 * min(step(xcoordinates(grid)), step(ycoordinates(grid)))
    ps = unique_points(ps, tol)

    Contour2D(ps, grid)
end


function _contour2D_grid(F::ContourFiber2D, grid::Grid2D,
    every_main_axis, samples_off_axis, options::MembershipTestOptions)

    ntries = options.ntries
    iterations = options.maxiters
    tol = options.tol

    fix_axis!(F, :w1)
    Fs = [deepcopy(F) for _=2:Threads.nthreads()]
    pushfirst!(Fs, F)

    B = Bitmap2D(grid)

    xmin, xmax, ymin, ymax = limits(grid)

    # FIRST PASS
    w1s = xcoordinates(grid)
    w2s = central_range(ymin, ymax, samples_off_axis)

    force_exit_y(x) = x[3] > ymax || x[3] < ymin

    Threads.@threads for k=1:every_main_axis:length(w1s)
        G = Fs[Threads.threadid()]
        w1 = w1s[k]
        update_fiber!(G, w1)
        px = xgridpoint(grid, w1)
        for w2 in w2s
            k = 0
            while k < ntries
                k += 1
                x = SVector(2π*rand(), 2π*rand(), w2)
                res = newton(G, x, iterations, tol, force_exit_y)
                if res.converged
                    py = ygridpoint(grid, res.solution[3])
                    safe_setindex!(B, true, py, px)
                end
            end
        end
    end

    # SECOND PASS
    for F in Fs
        fix_axis!(F, :w2)
    end
    w1s = central_range(xmin, xmax, samples_off_axis)
    w2s = ycoordinates(grid)

    force_exit_x(x) = x[3] > xmax || x[3] < ymin

    Threads.@threads for k=1:every_main_axis:length(w2s)
        G = Fs[Threads.threadid()]
        w2 = w2s[k]
        update_fiber!(G, w2)
        py = ygridpoint(grid, w2)
        for w1 in w1s
            k = 0
            while k < ntries
                k += 1
                x = SVector(2π*rand(), 2π*rand(), w1)
                res = newton(G, x, iterations, tol, force_exit_x)
                if res.converged
                    px = xgridpoint(grid, res.solution[3])
                    safe_setindex!(B, true, py, px)
                end
            end
        end
    end

    B
end


function isunique(vec::Vector{Float64}, x::Float64, tol::Float64)
    for y in vec
        if abs(x - y) < tol
            return false
        end
    end
    return true
end

# function unique_tol(vec::AbstractVector, tol::Float64)
#     out =
#     for y in vec
#         if abs(x - y) < tol
#             return false
#         end
#     end
#     return true
# end
