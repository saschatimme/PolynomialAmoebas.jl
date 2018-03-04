realified_coefficient_type(::Type{Complex{T}}) where T<:Real = promote_type(Float64, T)
realified_coefficient_type(::Type{T}) where T<:Real = promote_type(Float64, T)

"""
    realify(p::AbstractPolynomial)

Compute the realification of `p`. This is ``p(z) = p(x+iy) = p^{re}(x, y) + i p^{im}(x,y)``.
"""
function realify(p::MP.AbstractPolynomial, reals, imags)
    reals_imags = map((x, y) -> x + im*y, reals, imags)
    tmp = MP.subs(p, MP.variables(p) => reals_imags)

    parts = MP.map(MP.terms(tmp)) do term
        c = MP.coefficient(term)
        c_re = real(c)
        c_im = imag(c)
        (c_re * MP.monomial(term), c_im * MP.monomial(term))
    end

    f_re = MP.polynomial(first.(parts))
    f_im = MP.polynomial(last.(parts))

    f_re, f_im
end

function realify(p::MP.AbstractPolynomial)
    realify(p, realifyvariables(p)...)
end

function realifyvariables(p::MP.AbstractPolynomial)
    if MP.nvariables(p) == 2
        MP.@similarvariable p x1
        MP.@similarvariable p x2
        MP.@similarvariable p y1
        MP.@similarvariable p y2

        [x1, x2], [y1, y2]
    elseif MP.nvariables(p) == 3
        MP.@similarvariable p x1
        MP.@similarvariable p x2
        MP.@similarvariable p x3
        MP.@similarvariable p y1
        MP.@similarvariable p y2
        MP.@similarvariable p y3

        [x1, x2, x3], [y1, y2, y3]
    end
end
