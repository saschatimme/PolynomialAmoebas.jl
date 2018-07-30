"""
    AmoebaFiber2D(f::AbstractPolynomial, w=(0., 0.))

Construct the fiber `AmoebaFiber2D` of the bivariate polynomial `f` at `w`.
"""
function AmoebaFiber2D(p::MP.AbstractPolynomialLike{S}, w=(0.0, 0.0)) where S
    vars = MP.variables(p)
    @assert (length(vars) == 2) "Expected a bivariate polynomial, but got $(p)."

    T = realified_coefficient_type(S)
    f_re, f_im = SP.Polynomial.(1.0 .* realify(p))
    AmoebaFiber2D(f_re, f_im, w)
end

function AmoebaFiber2D(f_re::SP.Polynomial{T}, f_im::SP.Polynomial, w=(0.0, 0.0)) where T
    v = exp.(w)
    U = zeros(T, 2, 4)
    x = SVector{4, Float64}(0.0, 0.0, 0.0, 0.0)

    exponents = SP.exponents.((f_re, f_im))
    original_coefficients = copy.(SP.coefficients.((f_re, f_im)))

    AmoebaFiber2D(f_re, f_im, v, U, x, exponents, original_coefficients)
end

function update_fiber!(F::AmoebaFiber2D, w)
    w1, w2 = w
    F.v = (exp(w1), exp(w2))
    condition!(F)
    F
end


function condition!(fiber::AmoebaFiber2D)
    fiber.f_re.coefficients .= fiber.original_coefficients[1]
    fiber.f_im.coefficients .= fiber.original_coefficients[2]

    n1 = normalizer_amoeba_2d(fiber.exponents[1], fiber.f_re.coefficients, fiber.v)
    n2 = normalizer_amoeba_2d(fiber.exponents[2], fiber.f_im.coefficients, fiber.v)
    normalizer = inv(max(n1, n2))
    SP.scale_coefficients!(fiber.f_re, normalizer)
    SP.scale_coefficients!(fiber.f_im, normalizer)
    nothing
end

function normalizer_amoeba_2d(E, c, v)
    v1, v2 = v
    normalizer = sum(1:size(E, 2)) do j
        m1 = E[1,j]+E[3,j]
        m2 = E[2,j]+E[4,j]
        abs(c[j]) * Base.FastMath.pow_fast(v1, m1) *
        Base.FastMath.pow_fast(v2, m2)
    end
end


function phi!(fiber::AmoebaFiber2D, θ)
    v1, v2 = fiber.v
    x = SVector(v1 * cos(θ[1]), v2 * cos(θ[2]), v1 * sin(θ[1]), v2 * sin(θ[2]))
    fiber.x = x
    x
end

function evaluate(fiber::AmoebaFiber2D{T}, θ) where T
    x = phi!(fiber, θ)
    SVector(SP.evaluate(fiber.f_re, x), SP.evaluate(fiber.f_im, x))
end

function jacobian(fiber::AmoebaFiber2D, θ, precomputed=false)
    x = precomputed ? fiber.x : phi!(fiber, θ)
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
        x3, x4 = -x[3], -x[4]
        ∂re∂θ1 = U[1, 1] * x3
        ∂re∂θ2 = U[1, 2] * x4
        ∂im∂θ1 = U[2, 1] * x3
        ∂im∂θ2 = U[2, 2] * x4

        ∂re∂θ1 = muladd(U[1, 3], x[1], ∂re∂θ1)
        ∂re∂θ2 = muladd(U[1, 4], x[2], ∂re∂θ2)
        ∂im∂θ1 = muladd(U[2, 3], x[1], ∂im∂θ1)
        ∂im∂θ2 = muladd(U[2, 4], x[2], ∂im∂θ2)

        # SMatrix has column-major storage
        SMatrix{2,2}(∂re∂θ1, ∂im∂θ1, ∂re∂θ2, ∂im∂θ2)
    end

    out
end
