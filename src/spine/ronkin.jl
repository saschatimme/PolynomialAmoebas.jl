
"""
    ronkincoefficient(f, cc::ComponentComplement, N)

Compute the the Ronkin coefficient for the `ComponentComplement` `cc` using the trapezoidal rule
with `N` samples (which results in a steplength of ``2π/N``).
"""
function ronkincoefficient(f, cc::ComponentComplement, N)
    n_f = ronkinfunction(f, point(cc), N)

    n_f - order(cc)[1] * point(cc)[1] - order(cc)[2] * point(cc)[2]
end
function ronkincoefficients(f, ccs, nsamples=1024)
    map(cc -> ronkincoefficient(f, cc, nsamples), ccs)
end

"""
    ronkinfunction(f, z, N)

Evaluate the Ronkin function for `f` at `z` using the trapezoidal rule with `N` samples.

The Ronkin function is defined as
```math
x ⟼ \\frac{1}{(2πi)^2}∫_{\\text{Log}^{-1}|x|}\\log|f(z_1,z_2)| \\frac{dz_1dz_2}{z_1z_2}
```
where ``\\text{Log}^{-1}|x|`` is the fiber w.r.t. the log absolute map ``\\text{Log}|⋅|``.
Therefore
```math
\\text{Log}|(w_1,w_2)| = \\{ (exp(w_1+iθ_1), exp(w_2+iθ_2)) | θ_1, θ_2 ∈ [0,2π]\\}

Wit a change of variables ``(θ_1, ϴ_2) ⟼ (exp(w_1+iθ_1), exp(w_2+iθ_2))`` we get
```math
-\\frac{1}{(2πi)^2}∫_0^{2π}∫_0^{2π}\\log|f(exp(w_1+iθ_1), exp(w_2+iθ_2))|
    exp(w_1+iθ_1)exp(w_2+iθ_2) \\frac{dθ_1dθ_2}{exp(w_1+iθ_1)exp(w_2+iθ_2)}
```
which simplifies to
```math
\\frac{1}{2π}∫_0^{2π}\\frac{1}{2π}∫_0^{2π}\\log|f(exp(w_1+iθ_1), exp(w_2+iθ_2))| dθ_1dθ_2
```

Using the trapezoidal rule with `N` samples we can approximate the integral with
```math
1/N^2 ∑_{l=0}^{N-1}∑_{k=0}^{N-1} \\log|f(exp(w_1+iθ_1), exp(w_2+iθ_2))|
```
"""
function ronkinfunction(f, z, N)
    e_x, e_y = exp.(z)
    h = 2π/N
    # sum = 0.0
    partition = partition_work(N)
    inner_range = 0:N-1
    sums = Vector{Float64}(undef, Threads.nthreads())
    Threads.@threads for tid = 1:Threads.nthreads()
        r = shift_range(partition[tid], -1)
        partial_sum = ronkin_kernel(r, inner_range, f, h, e_x, e_y)
        sums[tid] = partial_sum
    end
    out = sum(sums) / N^2
    out
end

shift_range(r, n) = (r.start + n):(r.stop + n)

function ronkin_kernel(partial_outer_range, inner_range, f, h, e_x, e_y)
    partial_sum = 0.0
    for l in partial_outer_range, k in inner_range
        z1 = e_x * Base.FastMath.cis_fast(h*k)
        z2 = e_y * Base.FastMath.cis_fast(h*l)
        partial_sum += log(Base.FastMath.abs_fast(f(SVector(z1, z2))))
    end
    partial_sum
end

function ronkinpolynomial(p::MP.AbstractPolynomial, ronkincoeffs::Vector{Float64}, ccs)
    x, y = MP.variables(p)
    terms = map(ronkincoeffs, ccs) do r, cc
        i, j = order(cc)
        # We truncate the coefficients after 12 digits to avoid numerical problems
        r = trunc(r, digits=12)
        PolynomialAmoebas.Tropical(r) * x^i * y^j
    end
    MP.polynomial(terms)
end
