function mp_exponents_coefficients(poly::MP.AbstractPolynomialLike{T}, vars) where T
    terms = MP.terms(poly)
    nterms = length(terms)
    nvars = length(vars)
    exps = Vector{Vector{Int}}(nterms)
    coefficients = Vector{T}(nterms)
    for j = 1:nterms
        term = terms[j]
        coefficients[j] = MP.coefficient(term)
        exps[j] = [MP.degree(term, vars[i]) for i=1:nvars]
    end
    normalize_exponents_coeffs(exps, coefficients)
end

function normalize_exponents_coeffs(exponents::Matrix, coefficients::AbstractVector{T}) where T
    @assert length(coefficients) == size(exponents, 2)
    E = [exponents[:,j] for j=1:size(exponents, 2)]
    normalize_exponents_coeffs(E, coefficients)
end

function normalize_exponents_coeffs(E::Vector{Vector{Int}}, coefficients::AbstractVector{T}) where T
    P = sortperm(E, lt=revlexless, rev=true)
    f_exponents = Vector{Int}()
    f_coefficients = Vector{T}()
    k = 1
    while k ≤ length(P)
        v = E[P[k]]
        l = k
        while l < length(P) && E[P[l+1]] == v
            l += 1
        end
        append!(f_exponents, v)
        if l == k
            push!(f_coefficients, coefficients[P[k]])
        else
            c = sum(i -> coefficients[P[i]], k:l)
            push!(f_coefficients, c)

        end
        k = l + 1
    end
    f_exponents, f_coefficients
end

function revlexless(a, b)
    n = length(a)
    for j=n:-1:1
        if a[j] > b[j]
            return false
        elseif a[j] < b[j]
            return true
        end
    end
    return false
end


function convert_to_matrix(nvars, nterms, exponents)
    [exponents[nvars*(j - 1) + i] for i=1:nvars, j=1:nterms]
end

@inline pow(x::AbstractFloat, k::Integer) = Base.FastMath.pow_fast(x, k)
# simplified from Base.power_by_squaring
@inline function pow(x::Number, p::Integer)
    # if p == 1
    #     return copy(x)
    # elseif p == 0
    #     return one(x)
    # elseif p == 2
    #     return x*x
    # end
    t = trailing_zeros(p) + 1
    p >>= t
    while (t -= 1) > 0
        x *= x
    end
    y = x
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) >= 0
            x *= x
        end
        y *= x
    end
    return y
end

function static_pow(expr, k::Integer)
    if k == 2
        :($expr * $expr)
    elseif k == 3
        :(($expr * $expr * $expr))
    elseif k == 1
        :($expr)
    elseif k == 0
        :(one($expr))
    else
        :(pow($expr, $k))
    end
end


function x_(S::IntSet)
    str = ""
    for (i, a) in enumerate(S)
        if i == 1
            str *= "$(x_(a))"
        else
            str *= "_$(x_(a))"
        end
    end

    Symbol(str)
end
x_(i::Int) = Symbol("x", i)
x_(ik::NTuple{2, Int}) = ik[2] == 1 ? x_(ik[1]) : Symbol("x", ik[1], "_", ik[2])
u_(i::Int) = Symbol("u", i)
u_(i1::Int, i2::Int) = Symbol("u", i1, "_", i2)
u_(ik::NTuple{2, Int}) = ik[2] == 1 ? u_(ik[1]) : Symbol("u", ik[1], "_", ik[2])
c_(i::Int) = Symbol("c", i)

function occuring_exponents(M::Matrix, i::Int)
    exps = sort!(unique(M[i,:]))
    for j=2:length(exps)
        push!(exps, exps[j] - exps[j-1])
    end
    exps = sort!(unique(exps))
    if !isempty(exps) first(exps) == 0
        shift!(exps)
    end
    if !isempty(exps) && first(exps) == 1
        shift!(exps)
    end
    exps
end


"""
    group_powers(M::Matrix, NVars)

Create the most efficient way to compute all occuring powers.
This assumes that the input array is `x`. The new variables have the name
`xi_k` for ``x_i^k``.
"""
function group_powers(M::Matrix)
    xs = []
    for i = 1:size(M, 1)
        exps = occuring_exponents(M, i)
        push!(xs, :($(x_(i)) = x[$i]))
        if isempty(exps)
            continue
        end

        last_k = exps[1]
        last_xik = x_((i, last_k))
        p = static_pow(x_(i), last_k)
        push!(xs, :($last_xik = $p))
        for j=2:length(exps)
            k = exps[j]
            xik = x_((i,k))
            if _pow_already_computed(exps, k - last_k)
                push!(xs, :($xik = $last_xik * $(x_((i, k - last_k)))))
            else
                p = static_pow(x_(i), k - last_k)
                push!(xs, :($xik = $last_xik * $p))
            end
            last_k = k
            last_xik = xik
        end
    end
    return xs
end

function _pow_already_computed(exps, pow)
    for k in exps
        if k == pow
            return true
        elseif k > pow
            return false
        end
    end
    return false
end

"""
    batch_arithmetic_ops(op, operands)

Adds additional parenthes due to the following julia bug:
https://github.com/JuliaLang/julia/issues/14821
"""
function batch_arithmetic_ops(op::Symbol, ops)
    batch_size = 3
    l = 1
    if length(ops) == 1
        return :($(ops[1]))
    elseif length(ops) < batch_size + 2
        return Expr(:call, op, ops...)
    end
    batches = []
    while l ≤ length(ops)
        batch_end = min(length(ops), l + batch_size)
        if l == batch_end
            push!(batches, :($(ops[l])))
        else
            push!(batches, Expr(:call, op, ops[l:batch_end]...))
        end
        l = batch_end + 1
    end
    if length(batches) > 1
        return batch_arithmetic_ops(op, batches)
        # return Expr(:call, op, batches...)
    else
        return batches[1]
    end
end
