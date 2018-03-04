module SuperPolynomials
    import MultivariatePolynomials
    const MP = MultivariatePolynomials

    using Nullables

    export Polynomial

    struct StaticVal{T} end
    function Base.show(io::IO, ::Type{StaticVal{E}}) where {E}
        exps_hash = num2hex(hash(E))
        print(io, "StaticVal{$exps_hash}")
    end
    function Base.show(io::IO, ::StaticVal{E}) where {E}
        print(io, "StaticVal{$exps_hash}()")
    end


    """
        Polynomial([T, ] f::MP.AbstractPolynomial, [variables])
    """
    struct Polynomial{T, NVars, NTerms, Exponents<:StaticVal}
        coefficients::Vector{T}
    end

    function Polynomial(coefficients::Vector{T}, exponents::Vector{<:Integer}, NVars) where T
        @assert length(coefficients)*NVars == length(exponents)
        NTerms = length(coefficients)
        Exponents = StaticVal{ntuple(i -> exponents[i], length(exponents))}
        Polynomial{T, NVars, NTerms, Exponents}(coefficients)
    end

    function Polynomial(coefficients::Vector{T}, exponents::Matrix{<:Integer}) where T
        NVars = size(exponents, 1)
        exps, coeffs = normalize_exponents_coeffs(exponents, coefficients)
        Polynomial(coeffs, exps, NVars)
    end

    function Polynomial(p::MP.AbstractPolynomialLike, variables = MP.variables(p))
        exps, coeffs = mp_exponents_coefficients(p, variables)
        Polynomial(coeffs, exps, length(variables))
    end

    function Polynomial(::Type{T}, p::MP.AbstractPolynomialLike, variables = MP.variables(p)) where T
        exps, coeffs = mp_exponents_coefficients(p, variables)
        Polynomial(convert.(T, coeffs), exps, length(variables))
    end



    Base.show(io::IO, f::Polynomial) = show(io, typeof(f))


    include("superpolynomials/promotion_conversion.jl")
    include("superpolynomials/utilities.jl")
    include("superpolynomials/evaluate.jl")
    include("superpolynomials/horner.jl")
    include("superpolynomials/gradient_helpers.jl")
    include("superpolynomials/gradient.jl")

    export exponents, coefficients, scale_coefficients!

    """
        exponents(f)

    Get the exponents of `f` as a matrix where each column represents a monomial of the
    polynomial.
    """
    function exponents(::Polynomial{T, NVars, NTerms, StaticVal{Exponents}}) where {T, NVars, NTerms, Exponents}
        convert_to_matrix(NVars, NTerms, Exponents)
    end

    """
        coefficients(f)

    Get the coefficients of `f`.
    """
    coefficients(f::Polynomial) = f.coefficients

    """
        scale_coefficients!(f, λ)

    Scale the coefficients of `f` with the factor `λ`.
    """
    function scale_coefficients!(f::Polynomial, λ)
        scale!(f.coefficients, λ)
        f
    end
end
