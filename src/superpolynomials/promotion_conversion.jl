function Base.promote_rule(::Type{S}, ::Type{Polynomial{T, NV, NT, E}}) where {S<:Number, T, NV, NT, E}
    Polynomial{promote_type(S, T), NV, NT, E}
end
function Base.promote_rule(::Type{Polynomial{S, NV, NT, E}}, ::Type{Polynomial{T, NV, NT, E}}) where {S<:Number, T, NV, NT, E}
    Polynomial{promote_type(S, T), NV, NT, E}
end

function Base.convert(::Type{Polynomial{S, NV, NT, E}}, f::Polynomial{T, NV, NT, E}) where {S, T, NV, NT, E}
    Polynomial{S, NV, NT, E}(convert.(S, f.coefficients))
end
