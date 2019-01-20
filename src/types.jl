export AmoebaFiber2D, AmoebaFiber3D, CoamoebaFiber2D, CoamoebaFiber3D, ContourFiber2D,
    ImaginaryFiber2D, ImaginaryFiber3D

abstract type AbstractFiber{T, M, N} end

abstract type AbstractAmoebaFiber{T, N} <: AbstractFiber{T, 2, N} end
"""
An `AmoebaFiber2D` is a representation of the fiber of the amoeba ``\\mathcal{A}_f``
at ``(w_1,w_2)`` where `f` is a bivariate polynomial.
"""
mutable struct AmoebaFiber2D{T, E1, E2} <: AbstractAmoebaFiber{T, 2}
    f_re::SP.Polynomial{T, E1, Nothing}
    f_im::SP.Polynomial{T, E2, Nothing}

    # Tempory storage
    v::NTuple{2, Float64} # exp(w1), exp(w2)
    U::Matrix{T}
    x::SVector{4, Float64}

    # We use this to scale the system to improve the numerics
    exponents::NTuple{2, Matrix{Int}}
    original_coefficients::NTuple{2, Vector{T}}
end

"""
An `AmoebaFiber3D` is a representation of the fiber of the amoeba ``\\mathcal{A}_f``
at ``(w_1, w_2, w_3)`` where `f` is a trivariate polynomial.
"""
mutable struct AmoebaFiber3D{T, E1, E2} <: AbstractAmoebaFiber{T, 3}
    f_re::SP.Polynomial{T, E1, Nothing}
    f_im::SP.Polynomial{T, E2, Nothing}

    # Tempory storage
    v::NTuple{3, Float64} # exp(w1), exp(w2), exp(w3)
    U::Matrix{T}
    x::SVector{6, Float64}

    # We use this to scale the system to improve the numerics
    exponents::NTuple{2, Matrix{Int}}
    original_coefficients::NTuple{2, Vector{T}}
end

abstract type AbstractCoamoebaFiber{T, N} <: AbstractFiber{T, 2, N} end

"""
An `CoamoebaFiber2D` is a representation of the fiber of the coamoeba ``\\mathcal{A}^{\\prime}_f``
at ``(θ_1, θ_2)`` where `f` is a bivariate polynomial.
"""
mutable struct CoamoebaFiber2D{T, E1, E2} <: AbstractCoamoebaFiber{T, 2}
    f_re::SP.Polynomial{T, E1, Nothing}
    f_im::SP.Polynomial{T, E2, Nothing}

    # temporary storage
    sincosθ::NTuple{4, Float64}
    U::Matrix{T}
    x::SVector{4, Float64}
end

mutable struct CoamoebaFiber3D{T, E1, E2} <: AbstractCoamoebaFiber{T, 3}
    f_re::SP.Polynomial{T, E1, Nothing}
    f_im::SP.Polynomial{T, E2, Nothing}

    # temporary storage
    sincosθ::NTuple{6, Float64}
    U::Matrix{T}
    x::SVector{6, Float64}
end


abstract type AbstractContourFiber{T, N} <: AbstractFiber{T, 3, N} end
"""
An `ContourFiber2D` is a representation of the fiber of the contour of the amoeba ``\\mathcal{A}^{\\prime}_f``
at ``(w_1, w_2)`` along a onedimensional affine subspace
where `f` is a bivariate polynomial.
"""
mutable struct ContourFiber2D{T, E1, E2, E3} <: AbstractContourFiber{T, 3}
    f_re::SP.Polynomial{T, E1, Nothing}
    f_im::SP.Polynomial{T, E2, Nothing}
    g::SP.Polynomial{T, E3, Nothing}
    exponents::NTuple{3, Matrix{Int}}
    v_i::Float64
    x::SVector{4, Float64}
    U::Matrix{T}
    w1_fixed::Bool
    original_coefficients::NTuple{3, Vector{T}}
end


abstract type AbstractImaginaryFiber{T, N} <: AbstractFiber{T, 2, N} end

mutable struct ImaginaryFiber2D{T, E1, E2} <: AbstractImaginaryFiber{T, 2}
    f_re::SP.Polynomial{T, E1, Nothing}
    f_im::SP.Polynomial{T, E2, Nothing}

    y::NTuple{2, T}
    U::Matrix{T}
end

mutable struct ImaginaryFiber3D{T, E1, E2} <: AbstractImaginaryFiber{T, 3}
    f_re::SP.Polynomial{T, E1, Nothing}
    f_im::SP.Polynomial{T, E2, Nothing}

    y::NTuple{3, T}
    U::Matrix{T}
end


export AbstractAlgorithm, AbstractGridAlgorithm, ArchTrop, Simple, Greedy, Coarse

abstract type AbstractAlgorithm end
abstract type AbstractGridAlgorithm <: AbstractAlgorithm end

struct ArchTrop <: AbstractGridAlgorithm end
struct Coarse <: AbstractGridAlgorithm end
struct Simple <: AbstractGridAlgorithm end
struct Greedy <: AbstractGridAlgorithm end
struct Polygonal <: AbstractAlgorithm end


export NewtonPolygon
"""
    NewtonPolygon

Representation of a [Newton polygon](https://en.wikipedia.org/wiki/Newton_polygon) of a polynomial.

## Fields
- `polynomial::P`: the polynomial from which the NewtonPolygon was constructed.
- `lattices::Vector{SVector{2,Int}}`: the support of the polynomial `p`.
- `vertices::Vector{Int}`: the vertex indicdes of the Newton polygon.
- `facets::Vector{SVector{2,Int}}`: the facets of the Newton polygon, i.e. it's convex hull.
- `subdivision::Vector{SVector{3,Int}}`: the simplices of the regular subdivision of the
Newton polygon
"""
struct NewtonPolygon
    lattices::Vector{SVector{2, Int}}
    vertices::Vector{Int}
    facets::Vector{SVector{2, Int}}
    subdivision::Vector{Vector{Int}}
end


const Range = StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}

abstract type AbstractGrid{N} end
"""
A structure representing a 2D Grid.
"""
struct Grid2D <: AbstractGrid{2}
    xrange::Range
    yrange::Range
end

"""
A structure representing a 3D Grid.
"""
struct Grid3D <: AbstractGrid{3}
    xrange::Range
    yrange::Range
    zrange::Range
end

abstract type AbstractBitmap{N} end

export Bitmap2D, Bitmap3D
"""
A structure which holds a `BitMatrix` indicating active / inactive pixels and associated
world coordinates `xmin`, `xmax`, `ymin`, `ymax`.
"""
struct Bitmap2D <: AbstractBitmap{2}
    data::BitArray{2}
    grid::Grid2D
end

"""
Structute which holds a `BitArray{3}` indicating active / inactive voxels and associated
world coordinates `xmin`, `xmax`, `ymin`, `ymax`, `zmin`, `zmax`.
"""
struct Bitmap3D <: AbstractBitmap{3}
    data::BitArray{3}
    grid::Grid3D
end


export Tropical, TropicalCurve

"""
A `Tropical` number is an element of the semi-ring (ℝ ∪ {-∞}, ⊕, ⊙).
"""
struct Tropical{T<:Real} <: Number
    val::T
    isinf::Bool

    function Tropical{T}(val::T, isinf::Bool) where {T<:Real}
        new(val, isinf)
    end
    Tropical{T}(t::Tropical) where {T<:Real} = new(convert(T, t.val), t.isinf)
end

const AbstractTropicalPolynomial{T} = MP.AbstractPolynomial{Tropical{T}}

"""
A `TropicalCurve` is a tropical hypersurface of a bivariate polynomial.
"""
struct TropicalCurve{P<:AbstractTropicalPolynomial}
    polynomial::P
    vertices::Vector{SVector{2, Float64}}
    segments::Vector{SVector{2, Int64}}
    halfrays::Vector{Tuple{Int64, SVector{2, Float64}}}
end

export Spine2D, ComponentComplement, PolygonalAmoeba

"""
Represents a component of the complement of the amoeba.
It contains information about the order of the component, whether it is bounded
and one point out of the component.
"""
struct ComponentComplement{N}
    order::SVector{N, Int}
    x::SVector{N, Float64}
    bounded::Bool
end

"""
A `Spine2D` represents the spine of a two-dimensional amoeba. This also contains
informations about the components of the complement and the approximation used to obtain
this spine.
"""
struct Spine2D{P<:AbstractTropicalPolynomial}
    ronkin_polynomial::P
    curve::TropicalCurve{P}
    components_complement::Vector{ComponentComplement{2}}
    amoeba_approximation::Bitmap2D
end


"""
    PolygonalAmoeba

Holding the result of call to [amoeba](@ref) with the `Polygonal()` algorithm.
"""
struct PolygonalAmoeba{P}
    spine::Spine2D{P}
    polygons::Vector{NTuple{4, SVector{2, Float64}}}
    domain::NTuple{4, Float64}
    # qau
    error_bound::Float64
    iterations::Int
end

export Coamoeba
"""
    Coamoeba{N}

Holding a `BitArray` represeting the coamoeba embedded in ``[0,2π]^N``.
"""
struct Coamoeba{N}
    data::BitArray{N}
end


# We introduce some aliases
const Point = SVector{2, Float64}
const Direction = SVector{2, Float64}
const Pillar = Tuple{Point, Direction}
