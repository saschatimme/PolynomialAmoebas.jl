export ImaginaryFiber3D


"""
    ImaginaryFiber2D(f, y=(0., 0., 0.))

Construct the fiber `ImaginaryFiber3D` of the trivariate polynomial `f` at `y`.
"""
function ImaginaryFiber3D(p::MP.AbstractPolynomialLike{S}, y=(0.0, 0.0, 0.0)) where S
    vars = MP.variables(p)
    @assert (length(vars) == 3) "Expected a trivariate polynomial, but got $(p)."

    T = realified_coefficient_type(S)
    f_re, f_im = SP.Polynomial.( one(T) .* realify(p))

    ImaginaryFiber3D(f_re, f_im, float.(y), zeros(T, 2, 6))
end

function update_fiber!(F::ImaginaryFiber3D{T}, y) where T
    F.y = (y[1], y[2], y[3])
    F
end


function evaluate(fiber::ImaginaryFiber3D{T}, x::SVector{3,T}) where T
    xy = SVector(x[1], x[2], x[3], fiber.y[1], fiber.y[2], fiber.y[3])
    out = SVector(SP.evaluate(fiber.f_re, xy), SP.evaluate(fiber.f_im, xy))
    out
end

function jacobian(fiber::ImaginaryFiber3D, x::AbstractVector, precomputed=false)
    U = fiber.U
    xy = SVector(x[1], x[2], x[3], fiber.y[1], fiber.y[2], fiber.y[3])
    @inbounds begin
        ∇f_re = SP.gradient(fiber.f_re, xy)
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

    # SMatrix has column-major storage
    SMatrix{2, 3}(U[1,1], U[2,1], U[1, 2], U[2, 2], U[1, 3], U[2, 3])
end
