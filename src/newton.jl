struct NewtonResult{T, N}
    converged::Bool
    iterations::Int
    solution::SVector{N, T}
    residual::T
end

"""
    newton(F, x0, maxiterations, ftol)

Use a trust-region newton iteration[^1] to compute a root of the fiber `F` using maximal `iterations` iterations.
Convergence is declared if the infinity norm is less than `ftol`.

[^1]: Chapter 11.2 in Wright, Stephen J., and Jorge Nocedal. "Numerical optimization." Springer Science 35.67-68 (1999): 7.
"""
function newton(
    F::AbstractFiber{T, N, N},
    x0::SVector{N, T},
    iterations::Integer,
    ftol::Float64,
    abort=always_continue) where {N, T}

    factor = 1.0

    x = x0     # Current point

    r = evaluate(F, x)
    J = jacobian(F, x, true)

    it = 0
    norm_residual = maxnorm(r)
    converged = norm_residual < ftol

    Δ = factor * norm(x)
    if Δ == 0.0
        Δ = factor
    end

    η = 1e-4
    force_exit = false
    while !(converged || force_exit) && it < iterations
        it += 1
        # Compute proposed iteration step
        p, is_ok = dogleg(r, J, Δ)
        if !is_ok
            break
        end
        x += p

        r_new = evaluate(F, x)

        # Ratio of actual to predicted reduction (equation 11.45 in N&W)
        r_predict = r + J * p
        norm_r2 = sum(abs2,r)
        ρ = (norm_r2 - sum(abs2,r_new)) / (norm_r2 - sum(abs2,r_predict))

        if ρ > η
            # Successful iteration
            r = r_new
            J = jacobian(F, x, true)

            norm_residual = maxnorm(r)
            converged = norm_residual < ftol

            if abort(x)
                force_exit = true
            end
        else
            x -= p
            converged = false
        end

        # Update size of trust region
        if ρ < 0.1
            Δ *= 0.5
        elseif ρ >= 0.9
            Δ = 2 * norm(p)
        elseif ρ >= 0.5
            Δ = max(Δ, 2 * norm(p))
        end
    end

    NewtonResult(converged, it,  x, norm_residual)
end


always_continue(x) = false


"""
    dogleg(r, d, J, Δ::Real)
See Procedure 11.6 in "Numerical optimization."[^1] from Wright and Nocedal

[^1]: Wright, Stephen J., and Jorge Nocedal. "Numerical optimization." Springer Science 35.67-68 (1999): 7.
"""
function dogleg(r::SVector{N, T}, J::SMatrix{N, N, T}, Δ::Real) where {N,T}
    p_i = r
    # try
    p_i, is_ok = mysolve(J, r)
    if !is_ok
        return (p_i, false)
    end
    # catch
    #     # use Moore-Penrose inverse instead
    #     println("pinv")
    #     p_i = -SMatrix{N, N}(pinv(Matrix(J))) * r
    # end

    # Test if Gauss-Newton step is within the region
    if norm(p_i) <= Δ
        return p_i, true   # accepts equation 4.13 from N&W for this step
    else
        g = transposed_mul(J, r)

        # compute Cauchy point
        norm_g2 = sum(abs2, g)
        p_c = - norm_g2 / sum(abs2, J*g) .* g

        if norm(p_c) >= Δ
            # Cauchy point is out of the region, take the largest step along
            # gradient direction
            return g .* (-Δ / sqrt(norm_g2)), true
        else
            p_diff = p_i - p_c
            # Compute the optimal point on dogleg path
            b = 2 * dot(p_c, p_diff)
            a = dot(p_diff, p_diff)
            τ = (-b+sqrt(b^2 - 4a*(dot(p_c, p_c) - Δ^2)))/(2a)
            return p_c .+ τ .* p_diff, true
        end
    end
end

# We have a variant of Newton's method for underdetermined system
function newton(
    F::AbstractFiber{T, M, N},
    x0::SVector{N, T},
    iterations::Integer,
    ftol::Float64, abort=always_continue) where {M, N, T}

    x = x0     # Current point

    r = evaluate(F, x)
    norm_residual = maxnorm(r)
    J = jacobian(F, x, true)
    it = 0
    converged = norm_residual < ftol
    force_exit = false
    while !(converged || force_exit) && it < iterations
        it += 1
        # Compute proposed iteration step
        p, is_ok = mysolve(J, r)
        if !is_ok
            break
        end
        x -= p
        if abort(x)
            force_exit = true
        end

        r = evaluate(F, x)
        norm_residual = maxnorm(r)
        converged = norm_residual < ftol
    end

    NewtonResult(converged, it,  x, norm_residual)
end


function mysolve(M::SMatrix{2, 2, T}, u::SVector{2, T}) where T
    @inbounds minus_m3 = -M[3]
    @inbounds denom = muladd(M[1], M[4], minus_m3 * M[2])
    if abs(denom) < eps(T) * 10
        return u, false
    end
    invdenom = -inv(denom)
    @inbounds res = SVector(muladd(M[4], u[1], minus_m3 * u[2]) * invdenom, muladd(M[1], u[2], -M[2]*u[1]) * invdenom)
    return res, true
end

function mysolve(a::SMatrix{3, 3, T}, b::SVector{3, T}) where T
    d = det(a)
    if abs(d) < eps(T) * 10
        return b, false
        # throw(Base.LinAlg.SingularException(1))
    end
    invd = -inv(d)
    @inbounds out = SVector{3, T}(
        ((a[2,2]*a[3,3] - a[2,3]*a[3,2])*b[1] +
            (a[1,3]*a[3,2] - a[1,2]*a[3,3])*b[2] +
            (a[1,2]*a[2,3] - a[1,3]*a[2,2])*b[3]) * invd,
        ((a[2,3]*a[3,1] - a[2,1]*a[3,3])*b[1] +
            (a[1,1]*a[3,3] - a[1,3]*a[3,1])*b[2] +
            (a[1,3]*a[2,1] - a[1,1]*a[2,3])*b[3]) * invd,
        ((a[2,1]*a[3,2] - a[2,2]*a[3,1])*b[1] +
            (a[1,2]*a[3,1] - a[1,1]*a[3,2])*b[2] +
            (a[1,1]*a[2,2] - a[1,2]*a[2,1])*b[3]) * invd)
    return out, true
end

function mysolve(A::SMatrix{M, N, T}, u::SVector{M, T}) where {T, M, N}
    Q, R = qr(A')
    z = forward_substitution(R', u)
    Q * z, true
end

function forward_substitution(R::SMatrix{1, 1, T}, u::SVector{1, T}) where T
    SVector(u[1] / R[1, 1])
end
function forward_substitution(R::SMatrix{2, 2, T}, u::SVector{2, T}) where T
    x1 = u[1] / R[1, 1]
    x2 = (u[2] - R[2, 1] * x1) / R[2, 2]
    SVector(x1, x2)
end

function forward_substitution(R::SMatrix{3, 3, T}, u::SVector{3, T}) where T
    x1 = u[1] / R[1, 1]
    x2 = (u[2] - R[2, 1] * x1) / R[2, 2]
    x3 = (u[3] - R[3, 1] * x1 - R[3, 2] * x2) / R[3, 3]
    SVector(x1, x2, x3)
end

"""
    transposed_mul(J, r)

Compute `J'*r`.
"""
function transposed_mul(J::SMatrix{2, 2, T}, r::SVector{2, T}) where T
    @inbounds a = muladd(J[1, 1], r[1], J[2, 1] * r[2])
    @inbounds b = muladd(J[1, 2], r[1], J[2, 2] * r[2])
    SVector{2, T}(a, b)
end

function transposed_mul(J::SMatrix{3, 3, T}, r::SVector{3, T}) where T
    @inbounds a = J[1, 1] * r[1] + J[2, 1] * r[2] + J[3, 1] * r[3]
    @inbounds b = J[1, 2] * r[1] + J[2, 2] * r[2] + J[3, 2] * r[3]
    @inbounds c = J[1, 3] * r[1] + J[2, 3] * r[2] + J[3, 3] * r[3]
    SVector{3, T}(a, b, c)
end

function transposed_mul(J::SMatrix{2, 3, T}, r::SVector{2, T}) where T
    @inbounds a = J[1, 1] * r[1] + J[2, 1] * r[2]
    @inbounds b = J[1, 2] * r[1] + J[2, 2] * r[2]
    @inbounds c = J[1, 3] * r[1] + J[2, 3] * r[2]
    SVector{3, T}(a, b, c)
end


@inline maxnorm(x::SVector{2}) = max(abs(x[1]), abs(x[2]))
@inline maxnorm(x::SVector{3}) = max(abs(x[1]), abs(x[2]), abs(x[3]))
