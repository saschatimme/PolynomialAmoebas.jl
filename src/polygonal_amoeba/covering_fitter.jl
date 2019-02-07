const ReverseOrdering = Base.Order.ReverseOrdering{Base.Order.ForwardOrdering}

mutable struct CoveringFitter
    covering::Covering
    queue::PriorityQueue{Int, Float64, ReverseOrdering}
    approximated_global_error::Float64
    iteration::Int

    desired_maximal_error::Float64
    vertices_accuracy::Float64
    maximal_iterations::Int

    F::AmoebaFiber2D
    startvalue_generator::TorusStartValueGenerator{2}
    membership_test_options::MembershipTestOptions
end

function CoveringFitter(C::Covering, F::AmoebaFiber2D;
    desired_maximal_error=5e-2,
    vertices_accuracy = desired_maximal_error * 1e-3,
    maximal_iterations = 2000,
    options=MembershipTestOptions())

    covering = deepcopy(C)
    startvalue_generator = TorusStartValueGenerator{2}()
    initial_refinement!(covering, F, startvalue_generator, vertices_accuracy, options)

    queue = PriorityQueue{Int, Float64}(Base.Reverse)
    approximated_global_error = 0.0
    for (id, poly) in polygons(covering)
        loc_err = approximated_local_error(covering, poly, vertices_accuracy)
        approximated_global_error += loc_err
        enqueue!(queue, id => approximated_local_error(covering, poly, vertices_accuracy))
    end

    CoveringFitter(covering, queue, approximated_global_error, 0,
        desired_maximal_error, vertices_accuracy, maximal_iterations,
        F, startvalue_generator, options)
end

function update_error!(M::CoveringFitter, id::Int)
    if id != 0
        old_err = M.queue[id]
        new_err = approximated_local_error(M.covering, id, M.vertices_accuracy)
        M.queue[id] = new_err
        M.approximated_global_error -= (old_err - new_err)
    end
    nothing
end

function add2queue!(M::CoveringFitter, id::Int)
    if id != 0
        err = approximated_local_error(M.covering, id, M.vertices_accuracy)
        M.queue[id] = err
        M.approximated_global_error += err
    end
    nothing
end

"""
    fit!(covering_fitter)

Improve the covering until the `approximated_global_error` is less
than the `desired_maximal_error` or `maximal_iterations` were executed.
"""
function fit!(M::CoveringFitter; report_progress=false)
    k = 0
    while !done(M, k)
        k += 1
        if report_progress && mod(k, 5) == 0
            println("Iteration: $k - Error: $(M.approximated_global_error)")
        end
        next_step!(M)
    end
    k
end

function dequeue_pairs!(M::CoveringFitter, N::Int)
    ids = Int[]
    for k=1:N
        id, err = dequeue_pair!(M.queue)
        M.approximated_global_error -= err
        push!(ids, id)
    end
    ids
end

@inline function next_step!(M::CoveringFitter)
    M.iteration += 1
    id, err = dequeue_pair!(M.queue)
    M.approximated_global_error -= err
    pillar = central_pillar(M, id)
    poly = polygon(M.covering, id)

    id_left, id_right = split!(M.covering, id, pillar)
    update_error!(M, poly.left_id)
    update_error!(M, poly.right_id)
    add2queue!(M, id_left)
    add2queue!(M, id_right)
end

function Base.iterate(M::CoveringFitter, state=nothing)
    if state === nothing
        next_step!(M)
        M, state+1
    elseif done(M, state)
        return nothing
    else
        next_step!(M)
        M, state + 1
    end
end
done(M::CoveringFitter, k) = k > M.maximal_iterations || M.approximated_global_error < M.desired_maximal_error

Base.IteratorSize(::Type{<:CoveringFitter}) = Base.SizeUnknown()
Base.eltype(::Type{<:CoveringFitter}) = M

function initial_refinement!(covering::Covering, F::AbstractAmoebaFiber, startvalue_generator, ε, o::MembershipTestOptions)
    scalepillars!(covering) do origin, dir
        boundary_approximation(F, startvalue_generator, origin, dir, ε, o)
    end
    covering
end


"""
    boundary_approximation(F::AmoebaFiber, startvalue_generator, origin, direction, ε, membershiptest_options)

Starting from `origin` track the amoeba along `direction` with an initial step length.

Returns `t` s.t. `origin` + `t` * `direction` is *not* in the amoeba and the distance
A point is considered in the amoeba if we find a root (withi tolerance `membership_tol`) within `ntries` tries
and maximal `iterations` iterations. If we cannot find a root the step length is halfed.
"""
function boundary_approximation(F, startvalue_generator, origin::SVector{2}, direction::SVector{2}, ε::Float64, o::MembershipTestOptions, s0=0.0, Δs = 0.5)
    s = s0
    tol = ε / norm(direction)
    if tol ≥ 1.0
        return 1.0
    end
    while s ≤ 1.0
        w = origin + (s + Δs) * direction
        if membershiptest(F, w, startvalue_generator, o).successfull
            s += Δs
        else
            Δs *= 0.5
        end
        if Δs < tol
            return s + 2Δs
        end
    end
    return 1.0
end


function central_pillar(M::CoveringFitter, polygon_id::Int)
    poly = polygon(M.covering, polygon_id)
    orig, dir = average_pillar(M.covering, poly)
    pillar_line = (orig, orig + dir)
    l1, l2 = approximation_lines(M.covering, poly, M.vertices_accuracy)

    p1 = lineintersection(pillar_line, l1)
    p2 = lineintersection(pillar_line, l2)

    t1 = (p1 - orig)[1] / dir[1]
    t2 = (p2 - orig)[1] / dir[1]

    t1 = t1 > 1 ? 0.0 : t1
    t2 = t2 > 1 ? 0.0 : t2
    s0 = max(0.0, t1, t2)
    Δs = 0.5 * (1.0 - s0)

    t = boundary_approximation(M.F, M.startvalue_generator, orig, dir, M.vertices_accuracy, M.membership_test_options, s0, Δs)

    orig, t * dir
end

function average_pillar(covering::Covering, poly::CoveringPolygon)
    p1 = pillar(covering, poly.p1_id)
    p2 = pillar(covering, poly.p2_id)
    0.5 .* (p1 .+ p2)
end
