export ImaginaryFiber2D

"""
    ImaginaryFiber2D(f, y=(0., 0.))

Construct the fiber `ImaginaryFiber2D` of the bivariate polynomial `f` at `y`.
"""
function ImaginaryFiber2D(p::MP.AbstractPolynomialLike{S}, y=(0.0, 0.0)) where S
    vars = MP.variables(p)
    @assert (length(vars) == 2) "Expected a bivariate polynomial, but got $(p)."

    T = realified_coefficient_type(S)
    f_re, f_im = SP.Polynomial.( one(T) .* realify(p))

    ImaginaryFiber2D(f_re, f_im, y, zeros(T, 2, 4))
end

function update_fiber!(F::ImaginaryFiber2D{T}, y) where T
    F.y = (y[1], y[2])
    F
end


function evaluate(fiber::ImaginaryFiber2D{T}, x::SVector{2,T}) where T
    xy = SVector(x[1], x[2], fiber.y[1], fiber.y[2])
    out = SVector(SP.evaluate(fiber.f_re, xy), SP.evaluate(fiber.f_im, xy))
    out
end

function jacobian(fiber::ImaginaryFiber2D, x::AbstractVector, precomputed=false)
    U = fiber.U
    xy = SVector(x[1], x[2], fiber.y[1], fiber.y[2])
    @inbounds begin
        ∇f_re = SP.gradient(fiber.f_re, xy)
        for i=1:4
            U[1, i] = ∇f_re[i]
        end
        U[2, 3] = U[1, 1]
        U[2, 4] = U[1, 2]
        U[2, 1] = -U[1, 3]
        U[2, 2] = -U[1, 4]
    end

    # SMatrix has column-major storage
    SMatrix{2,2}(U[1,1], U[2,1], U[1, 2], U[2, 2])
end
