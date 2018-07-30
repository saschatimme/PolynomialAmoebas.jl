"""
    CoamoebaFiber2D(f::AbstractPolynomial, θ=(0., 0.))

Construct the fiber `CoamoebaFiber2D` of the bivariate polynomial `f` at `θ`.
"""
function CoamoebaFiber2D(p::MP.AbstractPolynomialLike{S}, θ=(0.0, 0.0)) where S
    vars = MP.variables(p)
    @assert (length(vars) == 2) "Expected a bivariate polynomial, but got $(p)."

    T = realified_coefficient_type(S)
    f_re, f_im = SP.Polynomial.(one(T) .* realify(p))

    sincosθ = _sincosθ(θ)
    U = zeros(2, 4)
    x = SVector(0.0, 0.0, 0.0, 0.0)
    CoamoebaFiber2D(f_re, f_im, sincosθ, U, x)
end

_sincosθ(θ) = (sin(θ[1]), sin(θ[2]), cos(θ[1]), cos(θ[2]))

function update_fiber!(F::CoamoebaFiber2D, θ)
    F.sincosθ = _sincosθ(θ)
    F
end

function phi!(fiber::CoamoebaFiber2D, w)
    s1, s2, c1, c2 = fiber.sincosθ
    v1, v2 = exp(w[1]), exp(w[2])
    x = SVector(v1 * c1, v2 * c2, v1 * s1, v2 * s2)
    fiber.x = x
    x
end

function evaluate(fiber::CoamoebaFiber2D, w)
    x = phi!(fiber, w)
    SVector(SP.evaluate(fiber.f_re, x), SP.evaluate(fiber.f_im, x))
end

function jacobian(fiber::CoamoebaFiber2D, w, precomputed=false)
    x = precomputed ? fiber.x : phi!(fiber, w)
    U = fiber.U
    @inbounds begin
        ∇f_re = SP.gradient(fiber.f_re, x)
        for i=1:4
            U[1, i] = ∇f_re[i]
        end
        U[2, 3] = U[1, 1]
        U[2, 4] = U[1, 2]
        U[2, 1] = -U[1, 3]
        U[2, 2] = -U[1, 4]
    end

    @inbounds out = begin
        ∂re∂w1 = muladd(U[1, 1], x[1], U[1, 3] * x[3])
        ∂re∂w2 = muladd(U[1, 2], x[2], U[1, 4] * x[4])
        ∂im∂w1 = muladd(U[2, 1], x[1], U[2, 3] * x[3])
        ∂im∂w2 = muladd(U[2, 2], x[2], U[2, 4] * x[4])

        SMatrix{2,2}(∂re∂w1, ∂im∂w1, ∂re∂w2, ∂im∂w2)
    end
    out
end
