"""
    ContourFiber2D(f; w1_fixed=true)

The fiber for the contour of `f` with the first axis fixed if `w1_fixed == true`.
"""
function ContourFiber2D(p::MP.AbstractPolynomialLike, w_i=0.0; w1_fixed=true)
    @assert (length(MP.variables(p)) == 2) "Expected a bivariate polynomial, but got $(p)."
    # MP.@similarvariable p s
    z1, z2 = MP.variables(p)
    realvars, imagvars = realifyvariables(p)

    fre, fim = realify(p, realvars, imagvars)
    g1re , g1im = realify(z1 * MP.differentiate(p, z1), realvars, imagvars)
    g2re , g2im = realify(z2 * MP.differentiate(p, z2), realvars, imagvars)
    g = g1re * g2im - g1im * g2re

    append!(realvars, imagvars)
    vars = realvars

    f_re = SP.Polynomial(1.0 * fre, vars)
    f_im = SP.Polynomial(1.0 * fim, vars)
    g = SP.Polynomial(1.0 * g, vars)

    v_i = exp(w_i)
    exponents = SP.exponents.((f_re, f_im, g))
    original_coefficients = copy.(SP.coefficients.((f_re, f_im, g)))
    ContourFiber2D(f_re, f_im, g, exponents, v_i, SVector(0.0, 0.0, 0.0, 0.0), zeros(3, 4), w1_fixed, original_coefficients)
end

function condition_contour_fiber!(f::SP.Polynomial, E, w1_fixed::Bool, v_i::Float64)
    if w1_fixed
        v1 = v_i
        m = maximum(map(j -> E[1,j]+E[3,j], 1:size(E, 2)))
        SP.scale_coefficients!(f, clamp(inv(v1^m), 1e-12, 1e12))
    else
        v2 = v_i
        m = maximum(map(j -> E[2,j]+E[4,j], 1:size(E, 2)))
        SP.scale_coefficients!(f, clamp(inv(v2^m), 1e-12, 1e12))
    end
end
function condition!(fiber::ContourFiber2D)
    fiber.f_re.coefficients .= fiber.original_coefficients[1]
    fiber.f_im.coefficients .= fiber.original_coefficients[2]
    fiber.g.coefficients .= fiber.original_coefficients[3]

    condition_contour_fiber!(fiber.f_re, fiber.exponents[1], fiber.w1_fixed, fiber.v_i)
    condition_contour_fiber!(fiber.f_im, fiber.exponents[2], fiber.w1_fixed, fiber.v_i)
    condition_contour_fiber!(fiber.g, fiber.exponents[3], fiber.w1_fixed, fiber.v_i)
end

function update_fiber!(c::ContourFiber2D, w_i::Float64)
    c.v_i = exp(w_i)
    condition!(c)
    c
end

function fix_axis!(c::ContourFiber2D, axis::Symbol)
    if axis == :w1
        c.w1_fixed = true
    elseif axis == :w2
        c.w1_fixed = false
    else
        throw(DomainError(":$axis is not a valid argument. Expected :w1 or :w2"))
    end
end

function phi!(fiber::ContourFiber2D, u)
    if fiber.w1_fixed
        v1 = fiber.v_i
        v2 = exp(u[3])
    else
        v1 = exp(u[3])
        v2 = fiber.v_i
    end
    x = SVector(v1 * cos(u[1]), v2 * cos(u[2]), v1 * sin(u[1]), v2 * sin(u[2]))
    fiber.x = x
    x
end

function evaluate(fiber::ContourFiber2D{T}, u::SVector{3,T}) where T
    x = phi!(fiber, u)

    SVector{3, T}(
        SP.evaluate(fiber.f_re, x),
        SP.evaluate(fiber.f_im, x),
        SP.evaluate(fiber.g, x))
end

function jacobian(fiber::ContourFiber2D, u::AbstractVector, precomputed=false)
    x = precomputed ? fiber.x : phi!(fiber, u)
    U = fiber.U
    @inbounds begin
        ∇f_re = SP.gradient(fiber.f_re, x)
        ∇f_im = SP.gradient(fiber.f_im, x)
        ∇g = SP.gradient(fiber.g, x)
        for i=1:4
            U[1, i] = ∇f_re[i]
            U[2, i] = ∇f_im[i]
            U[3, i] = ∇g[i]
        end
    end

    @inbounds out = begin
        mx3, mx4 = -x[3], -x[4]
        ∂f_re∂θ1 = muladd(U[1, 3], x[1], U[1, 1] * mx3)
        ∂f_re∂θ2 = muladd(U[1, 4], x[2], U[1, 2] * mx4)

        ∂f_im∂θ1 = muladd(U[2, 3], x[1], U[2, 1] * mx3)
        ∂f_im∂θ2 = muladd(U[2, 4], x[2], U[2, 2] * mx4)

        ∂g∂θ1 = muladd(U[3, 3], x[1], U[3, 1] * mx3)
        ∂g∂θ2 = muladd(U[3, 4], x[2], U[3, 2] * mx4)

        if fiber.w1_fixed
            ∂f_re∂w_i = muladd(U[1, 2], x[2], U[1, 4] * x[4])
            ∂f_im∂w_i = muladd(U[2, 2], x[2], U[2, 4] * x[4])
            ∂g∂w_i = muladd(U[3, 2], x[2], U[3, 4] * x[4])
        else
            ∂f_re∂w_i = muladd(U[1, 1], x[1], U[1, 3] * x[3])
            ∂f_im∂w_i = muladd(U[2, 1], x[1], U[2, 3] * x[3])
            ∂g∂w_i = muladd(U[3, 1], x[1], U[3, 4] * x[3])
        end

        # SMatrix has column-major storage
        SMatrix{3,3}(
            ∂f_re∂θ1, ∂f_im∂θ1, ∂g∂θ1,
            ∂f_re∂θ2, ∂f_im∂θ2, ∂g∂θ2,
            ∂f_re∂w_i, ∂f_im∂w_i, ∂g∂w_i)
    end
    out
end
