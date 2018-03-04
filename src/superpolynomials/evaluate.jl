export evaluate

function evaluate_impl(::Type{T}, M::Matrix{Int}) where T
    m, n = size(M)
    grouped_powers = group_powers(M)

    as = [:(out = zero($T))]
    for j = 1:n
        factors = []
        push!(factors, :(f.coefficients[$j]))
        nmultiplications = sum(M[:, j] .> 0)
        mult_counter = 0
        if nmultiplications == 0
            push!(as, :(@inbounds out += $(factors[1])))
        else
            for i=1:m
                k = M[i, j]
                if k > 0
                    mult_counter += 1
                    xik = x_((i, k))
                    if mult_counter == nmultiplications
                        a = batch_arithmetic_ops(:*, factors)
                        push!(as, :(@inbounds out = muladd($a, $xik, out)))
                    else
                        push!(factors, :($xik))
                    end
                end
            end
        end
    end

    Expr(:block,
        grouped_powers...,
        as...,
        :(out))
end

function evaluate_impl(::Type{T}, f::Type{Polynomial{S, NVars, NTerms, StaticVal{Exponents}}}) where {T, S, NVars, NTerms, Exponents}
    M = convert_to_matrix(NVars, NTerms, Exponents)
    evaluate_impl(promote_type(T, S), M)
end

@generated function evaluate(f::Polynomial{S, NVars, NTerms, StaticVal{Exponents}}, x::AbstractVector{T}) where {T, S, NVars, NTerms, Exponents}
    evaluate_impl(T, f)
end

@inline (f::Polynomial)(x) = evaluate(f, x)
