export inf, tropicalpolynomial

import Base: +, *, ==, ^, -

"""
    Tropical(val [, isinf=false])

Create a number in the Tropical semi ring where ``a ⨁ b = max(a,b)`` and ``a ⨀ b = a + b``.
"""
Tropical(val::T, isinf=false) where {T<:Real} = Tropical{T}(val, isinf)
inf(val::T) where {T<:Real} = Tropical(val, true)
inf(::Type{T}) where {T<:Real} = Tropical(zero(T), true)
inf(t::Tropical{T}) where {T<:Real} = Tropical(t.val, true)

Base.isinf(t::Tropical) = t.isinf

function (+)(a::Tropical{T}, b::Tropical{T}) where {T<:Real}
    if isinf(a)
        b
    elseif isinf(b)
        a
    else
        Tropical(max(a.val, b.val))
    end
end

# ⨁(a::Tropical{T}, b::Tropical{T}) where {T<:Real} = a + b
# ⨁(a::Tropical, b::Tropical) = +(promote(a,b)...)
# ⨁(a::T, b::T) where {T<:Real} = Tropical(max(a, b))
# ⨁(a::Real, b::Real) = ⨁(promote(a, b)...)

# this doesn't make sense, but we define it nonetheless to make MultivariatePolynomials
# happy
(-)(a::Tropical) = a

function (*)(a::Tropical{T}, b::Tropical{T}) where {T<:Real}
    if isinf(a)
        a
    elseif isinf(b)
        b
    else
        Tropical(a.val + b.val)
    end
end

# ⨀(a::Tropical{T}, b::Tropical{T}) where {T<:Real} = a * b
# ⨀(a::Tropical, b::Tropical) = *(promote(a,b)...)
# ⨀(a::T, b::T) where {T<:Real} = Tropical(a + b)
# ⨀(a::Real, b::Real) = ⨀(promote(a, b)...)


(^)(a::Tropical, p::Integer) = isinf(a) ? a : Tropical(p * a.val)

function Base.isequal(a::Tropical{T}, b::Tropical{T}) where T
    if isinf(a) || isinf(b)
        isinf(a) == isinf(b)
    else
        a.val == b.val
    end
end
(==)(a::Tropical{T}, b::Tropical{T}) where T = isequal(a,b)


Base.one(::Tropical{T}) where T = Tropical(zero(T))
Base.one(::Type{Tropical{T}}) where T = Tropical(zero(T))

Base.zero(::Tropical{T}) where {T} = Tropical(one(T), true)
Base.zero(::Type{Tropical{T}}) where {T} = Tropical(one(T), true)

Base.show(io::IO, t::Tropical{T}) where {T} = isinf(t) ? print(io, "-∞") : print(io, t.val)

function Base.promote_rule(::Type{Tropical{T}}, ::Type{Tropical{S}}) where {S,T}
    Tropical{promote_type(S,T)}
end
function Base.promote_rule(::Type{T}, ::Type{Tropical{S}}) where {S<:Real,T<:Real}
    Tropical{promote_type(S,T)}
end



Base.convert(::Type{Tropical{T}}, t::Tropical{T}) where {T<:Real} = t
function Base.convert(::Type{Tropical{T}}, t::Tropical{<:Real}) where {T<:Real}
    Tropical{T}(convert(T, t.val), t.isinf)
end
function Base.convert(::Type{Tropical{T}}, x::Real) where {T<:Real}
    if _isone(x)
        one(Tropical{T})
    elseif iszero(x)
        zero(Tropical{T})
    else
        Tropical(convert(T, x))
    end
end

# function Base.convert(::Type{T}, t::Tropical{<:Real}) where {T<:Number}
#     if _isone(t)
#         one(T)
#     elseif iszero(t)
#         zero(T)
#     else
#         convert(T, t.val)
#     end
# end


"""
    tropicalpolynomial(p::MP.AbstractPolynomial)

Interprets a normal polynomials a tropical one. This will especially convert all coefficents which are
1 (neutral element w.r.t. multiplication) to 0 (neutral element wr.t. tropical multiplication).
"""
function tropicalpolynomial(p::MP.AbstractPolynomial)
    mapcoefficients(c -> convert(Tropical{Float64}, c), p)
end


Base.float(t::Tropical) = isinf(t) ? -Inf : float(t.val)
