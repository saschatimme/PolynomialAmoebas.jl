export MembershipTestResult, MembershipTestOptions, membershiptest

"""
    MembershipTestResult{T, N}

**Fields**
* `successfull::Bool`
* `tries::Int` The number of times Newton's method was started
* `converged_iterations::Int` Number iterations needed in a successfull iteration. If the test was not successfull `imax` is the maximal number of iterations.
* `startvalue::SVector{N, T}` The start value of the successfull iteration. Otherwise a random value.
* `solution::SVector{N, T}` The solution value of the successfull iteration. Otherwise a random value.
"""
struct MembershipTestResult{T, N}
    successfull::Bool
    tries::Int
    converged_iterations::Int
    startvalue::SVector{N, T}
    solution::SVector{N, T}
end


"""
    MembershipTestOptions(iterations=50, ntries=30, tol=1e-7)

# Details
* `iterations` The maximal number of iterations in one try.
* `ntries` The number of times Newton's method will be tried
* `tol` The desired accuracy when Newton's method converges.
"""
struct MembershipTestOptions
    maxiters::Int
    tol::Float64
    ntries::Int
end
function MembershipTestOptions(;maxiters=50, ntries=30, tol=1e-7)
    MembershipTestOptions(maxiters, tol, ntries)
end


"""
    membershiptest(F::AbstractFiber{T}, w:, options=MembershipTestOptions()) where T

Check wether the point `w` is contained in the fiber `F`.
The tests uses Newton's method with configuration provided by `options`.
Returns a [`MembershipTestResult`](@ref).
"""
function membershiptest(F, w, gen, o=MembershipTestOptions())
    membershiptest(F, w, gen, (startvalue(gen),), o)
end
function membershiptest(F::AbstractFiber{T}, w, gen::StartValueGenerator, initial_startvalues, o::MembershipTestOptions) where T
    update_fiber!(F, w)
    k = 0
    i = 1
    if isempty(initial_startvalues)
        θ = startvalue(gen)
    else
        θ = initial_startvalues[i]
    end

    maxiters = o.maxiters
    while true
        k += 1
        # try
        res = newton(F, θ, maxiters, o.tol)
        if res.converged
            return MembershipTestResult(true, k, res.iterations, θ, res.solution)
        # If it seems that we could get a zero we iterate a little bit more
        # We should only need a small amount of iterations due to the
        # locally quadratically convergence of the newton iteration
        elseif res.residual < 1e-2
            refined = newton(F, res.solution, 10, o.tol)
            if refined.converged
                return MembershipTestResult(true, k, maxiters + refined.iterations, θ, res.solution)
            end
        end
        # catch err
        #     @show err
        # end
        if k > o.ntries
            break
        end
        if i < length(initial_startvalues)
            i += 1
            θ = initial_startvalues[i]
        else
            maxiters = o.maxiters
            θ = startvalue(gen)
        end
    end
    return MembershipTestResult(false, o.ntries, maxiters, θ, θ)
end
