"""
    CoamoebaFiber3D(f::AbstractPolynomial, θ=(0., 0., 0.))

Construct the fiber `CoamoebaFiber3D` of the trivariate polynomial `f` at `θ`.
"""
function CoamoebaFiber3D(p::MP.AbstractPolynomialLike{S}, θ=(0.0, 0.0, 0.0)) where S
    vars = MP.variables(p)
    @assert (length(vars) == 3) "Expected a trivariate polynomial, but got $(p)."

    T = realified_coefficient_type(S)
    f_re, f_im = SP.Polynomial.(one(T) .* realify(p))

    sincosθ = _sincosθ3(θ)
    U = zeros(2, 6)
    x = SVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    CoamoebaFiber3D(f_re, f_im, sincosθ, U, x)
end

_sincosθ3(θ) = (sin(θ[1]), sin(θ[2]), sin(θ[3]), cos(θ[1]), cos(θ[2]), cos(θ[3]))

function update_fiber!(F::CoamoebaFiber3D, θ)
    F.sincosθ = _sincosθ3(θ)
    F
end

function phi!(fiber::CoamoebaFiber3D, w)
    s1, s2, s3, c1, c2, c3 = fiber.sincosθ
    v1, v2, v3 = exp.(w)
    x = SVector(v1 * c1, v2 * c2, v3 * c3, v1 * s1, v2 * s2, v3 * s3)
    fiber.x = x
    x
end

function evaluate(fiber::CoamoebaFiber3D, w)
    x = phi!(fiber, w)
    SVector(SP.evaluate(fiber.f_re, x), SP.evaluate(fiber.f_im, x))
end

function jacobian(fiber::CoamoebaFiber3D, w, precomputed=false)
    x = precomputed ? fiber.x : phi!(fiber, w)
    U = fiber.U
    @inbounds begin
        ∇f_re = SP.gradient(fiber.f_re, x)
        for i=1:6
            U[1, i] = ∇f_re[i]
        end
        U[2, 4] = U[1, 1]
        U[2, 5] = U[1, 2]
        U[2, 6] = U[1, 3]
        U[2, 1] = -U[1, 4]
        U[2, 2] = -U[1, 5]
        U[2, 3] = -U[1, 6]
    end

    @inbounds out = begin
        ∂re∂w1 = muladd(U[1, 1], x[1], U[1, 4] * x[4])
        ∂re∂w2 = muladd(U[1, 2], x[2], U[1, 5] * x[5])
        ∂re∂w3 = muladd(U[1, 3], x[3], U[1, 6] * x[6])
        ∂im∂w1 = muladd(U[2, 1], x[1], U[2, 4] * x[4])
        ∂im∂w2 = muladd(U[2, 2], x[2], U[2, 5] * x[5])
        ∂im∂w3 = muladd(U[2, 3], x[3], U[2, 6] * x[6])

        SMatrix{2, 3}(∂re∂w1, ∂im∂w1, ∂re∂w2, ∂im∂w2, ∂re∂w3, ∂im∂w3)
    end
    out
end
