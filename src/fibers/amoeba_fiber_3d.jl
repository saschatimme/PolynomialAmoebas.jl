"""
    AmoebaFiber3D(f::AbstractPolynomial, w=(0., 0., 0.))

Construct the fiber `AmoebaFiber3D` of the bivariate polynomial `f` at `w`.
"""
function AmoebaFiber3D(p::MP.AbstractPolynomialLike{S}, w=(0.0, 0.0, 0.0)) where S
    vars = MP.variables(p)
    @assert (length(vars) == 3) "Expected a trivariate polynomial, but got $(p)."

    T = realified_coefficient_type(S)
    f_re, f_im = SP.Polynomial.(1.0 .* realify(p))

    v = exp.(w)
    U = zeros(T, 2, 6)
    x = SVector{6, Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    exponents = SP.exponents.((f_re, f_im))
    original_coefficients = copy.(SP.coefficients.((f_re, f_im)))

    AmoebaFiber3D(f_re, f_im, v, U, x, exponents, original_coefficients)
end

function update_fiber!(F::AmoebaFiber3D, w)
    w1, w2, w3 = w
    F.v = (exp(w1), exp(w2), exp(w3))
    condition!(F)
    F
end

function condition!(fiber::AmoebaFiber3D)
    fiber.f_re.coefficients .= fiber.original_coefficients[1]
    fiber.f_im.coefficients .= fiber.original_coefficients[2]

    n1 = normalizer_amoeba_3d(fiber.exponents[1], fiber.v)
    n2 = normalizer_amoeba_3d(fiber.exponents[2], fiber.v)
    normalizer = inv(max(n1, n2))
    SP.scale_coefficients!(fiber.f_re, normalizer)
    SP.scale_coefficients!(fiber.f_im, normalizer)
    nothing
end

function normalizer_amoeba_3d(E, v)
    v1, v2, v3 = v
    normalizer = sum(1:size(E, 2)) do j
        m1 = E[1,j]+E[4,j]
        m2 = E[2,j]+E[5,j]
        m3 = E[3,j]+E[6,j]
        Base.FastMath.pow_fast(v1, m1) *
        Base.FastMath.pow_fast(v2, m2) *
        Base.FastMath.pow_fast(v3, m3)
    end
end

function condition_amoeba_fiber_3d!(f, E, v)
    v1, v2, v3 = v
    normalizer = sum(1:size(E, 2)) do j
        m1 = E[1,j]+E[4,j]
        m2 = E[2,j]+E[5,j]
        m3 = E[3,j]+E[6,j]
        Base.FastMath.pow_fast(v1, m1) *
        Base.FastMath.pow_fast(v2, m2) *
        Base.FastMath.pow_fast(v3, m3)
    end
    SP.scale_coefficients!(f, 1 / normalizer)
end


function phi!(fiber::AmoebaFiber3D, θ)
    v1, v2, v3 = fiber.v
    x = SVector(
        v1 * cos(θ[1]), v2 * cos(θ[2]), v3 * cos(θ[3]),
        v1 * sin(θ[1]), v2 * sin(θ[2]), v3 * sin(θ[3])
        )
    fiber.x = x
    x
end

function evaluate(fiber::AmoebaFiber3D, θ)
    x = phi!(fiber, θ)
    SVector(SP.evaluate(fiber.f_re, x), SP.evaluate(fiber.f_im, x))
end

function jacobian(fiber::AmoebaFiber3D, θ, precomputed=false)
    x = precomputed ? fiber.x : phi!(fiber, θ)
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
        x4, x5, x6 = -x[4], -x[5], -x[6]
        ∂re∂θ1 = U[1, 1] * x4
        ∂re∂θ2 = U[1, 2] * x5
        ∂re∂θ3 = U[1, 3] * x6
        ∂im∂θ1 = U[2, 1] * x4
        ∂im∂θ2 = U[2, 2] * x5
        ∂im∂θ3 = U[2, 3] * x6

        ∂re∂θ1 = muladd(U[1, 4], x[1], ∂re∂θ1)
        ∂re∂θ2 = muladd(U[1, 5], x[2], ∂re∂θ2)
        ∂re∂θ3 = muladd(U[1, 6], x[3], ∂re∂θ3)
        ∂im∂θ1 = muladd(U[2, 4], x[1], ∂im∂θ1)
        ∂im∂θ2 = muladd(U[2, 5], x[2], ∂im∂θ2)
        ∂im∂θ3 = muladd(U[2, 6], x[3], ∂im∂θ3)

        # SMatrix has column-major storage
        SMatrix{2,3}(∂re∂θ1, ∂im∂θ1, ∂re∂θ2, ∂im∂θ2, ∂re∂θ3, ∂im∂θ3)
    end

    out
end
