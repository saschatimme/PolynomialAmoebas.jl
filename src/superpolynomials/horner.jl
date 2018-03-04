export evalhorner

function next_k(I, n, d_I)
    k = 1
    for j = d_I:-1:1
        if I[n][j] != I[n-1][j]
            return j
        end
    end
    k
end

function decrease_lex(a)
    b = copy(a)
    for i=1:length(a)
        if a[i] > 0
            b[i] -= 1
            return b
        end
    end
    b
end


function lower_set(E)
    sorted = sort(E, lt=revlexless, rev=true)
    zero_el = zero(last(sorted))
    if zero_el != last(sorted)
        push!(sorted, zero_el)
    end
    set = [shift!(sorted)]
    while !isempty(sorted)
        next = shift!(sorted)
        next_lower = decrease_lex(last(set))
        while next_lower != next && !revlexless(next_lower, next)
            push!(set, next_lower)
            next_lower = decrease_lex(next_lower)
        end
        push!(set, next)
    end
    set
end


function horner_impl(f::Type{Polynomial{T, NVars, NTerms, StaticVal{Exponents}}}) where {T, NVars, NTerms, Exponents}
    matrix = convert_to_matrix(NVars, NTerms, Exponents)
    E = [matrix[:,j] for j=1:NTerms]
    P = sortperm(E, lt=revlexless, rev=true)
    I = lower_set(E)

    d_I = NVars
    M = length(I)

    # we save all powers computed to precompute them later
    xpowers = [Vector{Int}() for _=1:d_I]

    _xs = []
    push!(_xs, :(coefficients = f.coefficients))
    push!(_xs, :(@inbounds r0 = coefficients[$(P[1])]))
    r0_zero = false
    last_k = 0
    l = 2
    r = [Symbol("r", k) for k=1:d_I]
    rs_zero = trues(d_I)

    for k=1:d_I
        push!(_xs, :($(r[k]) = zero($T)))
    end

    n = 2
    while n ≤ M
        # k = max{1≤j≤d_I : α(n)_j ≠ α(n-1)_j}
        k = next_k(I, n, d_I)
        n_i = n + 1

        # We want to find out if we can write instead of
        #
        # r1 = x1 * (r0 + r1)
        # r0 = zero(Float64)
        # r1 = x1 * (r0 + r1)
        #
        # r1 = x1_2 * (r0 + r1)
        if l ≤ NTerms && I[n] != E[P[l]]
            while n_i ≤ M
                n_k = next_k(I, n_i, d_I)
                # @show n_k
                if n_k != k
                    break
                end
                if l ≤ NTerms && (I[n_i] == E[P[l]])
                    n_i += 1
                    break
                end
                n_i += 1
            end
        end

        non_zero_r = r[filter(i -> !rs_zero[i], 1:k)]
        if r0_zero && length(non_zero_r) > 1
            summand = batch_arithmetic_ops(:+, non_zero_r)
            # summand = Expr(:call, :+, non_zero_r...)
        elseif r0_zero && isempty(non_zero_r)
            summand = :(zero($T))
        elseif r0_zero
            summand = :($(non_zero_r[1]))
        elseif isempty(non_zero_r)
            summand = :r0
        else
            e = batch_arithmetic_ops(:+, non_zero_r)
            summand = :(r0 + ($e))
            # summand = Expr(:call, :+, :r0, non_zero_r...)
        end
        if n_i - n > 1
            push!(xpowers[k], n_i - n)
            xkj = Symbol("x", k, "_", (n_i - n))
            push!(_xs, :($(r[k]) = $(xkj)*($(summand))))
            n = n_i - 1
        else
            xk = Symbol("x", k)
            push!(_xs, :($(r[k]) = $xk*($(summand))))
        end
        rs_zero[k] = false
        if l ≤ NTerms && I[n] == E[P[l]]
            push!(_xs, :(@inbounds r0 = coefficients[$(P[l])]))
            r0_zero = false
            l += 1
        else
            if n == M
                push!(_xs, :(r0 = zero($T)))
            end
            r0_zero = true
        end
        for i = 1:k-1
            if !rs_zero[i] && n != M
                push!(_xs, :($(r[i]) = zero($T)))
            end
            rs_zero[i] = true
        end
        n += 1
    end

    # We now still have to compute the powers
    _ps = []
    for (i, x_i_powers) in enumerate(xpowers)
        i_exps = sort!(unique(x_i_powers))
        if !isempty(i_exps)
            last_k = first(i_exps)
            xi = Symbol("x", i)
            push!(_ps, :($xi = x[$i]))
            last_xik = Symbol("x", i, "_", last_k)
            if last_k > 1
                push!(_ps, :($last_xik = $(static_pow(xi, last_k))))
            else
                push!(_ps, :($last_xik = $xi))
            end
            for k in Iterators.drop(i_exps, 1)
                xik = Symbol("x", i, "_",k)
                if (k - last_k) > 1
                    push!(_ps, :($xik = $last_xik * pow($xi, $(k - last_k))))
                else
                    push!(_ps, :($xik = $last_xik * $xi))
                end
                last_k = k
                last_xik = xik
            end
        else
            xi = Symbol("x", i)
            push!(_ps, :($xi = x[$i]))
        end
    end

    e = batch_arithmetic_ops(:+, r[.!rs_zero])
    out = :(r0 + ($e))
    Expr(:block,
        _ps...,
        _xs...,
        out)
end

"""
    evalhorner(f, x)

Evaluate `f(x)` using a multivariate extension to Horner's methods.
For some background and error estimates see "On the multivariate Horner scheme." from
J. Peña[^1].

[^1]: Peña, Juan Manuel. "On the multivariate Horner scheme." SIAM journal on numerical analysis 37.4 (2000): 1186-1197.
"""
@generated function evalhorner(f::Polynomial{T, NVars, NTerms, StaticVal{Exponents}}, x) where {T, NVars, NTerms, Exponents}
    horner_impl(f)
end
