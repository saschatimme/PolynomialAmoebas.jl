export spine, components_complement, ronkin_polynomial, hypersurface, amoeba_approximation

function Spine2D(polynomial::AbstractTropicalPolynomial, ccs::Vector{<:ComponentComplement}, aa::Bitmap2D)
    Spine2D(polynomial, TropicalCurve(polynomial), ccs, aa)
end

Spine2D(f::MP.AbstractPolynomial; kwargs...) = spine(f; kwargs...)

"""
    spine(f::MP.AbstractPolynomial; options...)

Compute the spine of the amoeba ``\\mathcal{A}(f)``. This algorithm computes first an approximation of
``\\mathcal{A}(f)`` and from this the spine. Returns a Spine2D.

## Example
```julia
@polyvar x y
# use all the defaults
spine(x^2 + y^2 + 1)

# Maybe we think that the compleement of the amoeba has some very small components
# and we want to use an explicit domain
spine(x^2 + y^2 + 1, domain=(-5, 5, -5, 5), minimal_component_size=0.001)
```

Optional arguments:
* `minimal_component_size=0.01`: A guarantee that in each component of the complement of ``\\mathcal{A}(f)`` fits a ball with this diameter.
If this does not hold the algorithm can return a wrong result.
* `domain`: A tuple in the form `(xmin, xmax, ymin, ymax)` which defines a section Ω
for which the amoeba ``\\mathcal{A}(f)`` is computed. This domain has to be such that the intersection
Ω ∩ ``\\mathcal{A}(f)`` still captures the correct topology of ``\\mathcal{A}(f)``.
* `grid`: Based on `minimal_component_size` and `domain` a grid can be computed automatically.
This can also be overwritten with this option.
* `membership_options=[MembershipTestOptions](@ref)()`: Options for the membership test.
* `nsamples=1024`: For the computation we numerically evaluate intergrals. This is the number
of sample points used in the computation.


    spine(f::MP.AbstractPolynomial, A::Bitmap2D; nsamples=1024)

Compute the spine based on the given amoeba approximation `A`.
"""
function spine(f::MP.AbstractPolynomial;
    domain=amoeba_carcase_domain_heuristic(f),
    minimal_component_size=0.01,
    grid = spine_grid(domain..., minimal_component_size),
    membership_options=MembershipTestOptions(),
    kwargs...)
    A = amoeba(f, alg=Greedy(), grid=grid, membership_options=membership_options)
    spine(f, A; kwargs...)
end

function spine(f::MP.AbstractPolynomial, A::Bitmap2D; nsamples=1024)
    p = SP.Polynomial((1.0+0.0im)*f)
    cc, _ = components_complement(A, AmoebaFiber2D(f), p, nsamples=nsamples)

    ronkincoeffs = ronkincoefficients(p, cc)
    ronkinpoly = ronkinpolynomial(f, ronkincoeffs, cc)
    Spine2D(ronkinpoly, cc, A)
end

function spine_grid(f, δ; kwargs...)
    spine_grid(amoeba_carcase_domain_heuristic(f; kwargs...)..., δ)
end
function spine_grid(xmin, xmax, ymin, ymax, δ)
    rx = ceil(Int, (xmax-xmin) / √(0.5δ))
    ry = ceil(Int, (ymax-ymin) / √(0.5δ))
    Grid2D(xmin, xmax, ymin, ymax, rx, ry)
end

"""
    components_complement(S::Spine2D)

Get the ComponentComplement of the spine `S`.
"""
components_complement(S::Spine2D) = S.components_complement


"""
    ronkin_polynomial(S::Spine2D)

Get the tropical polynomial which defines this hypersurface.
"""
ronkin_polynomial(S::Spine2D) = S.ronkin_polynomial


"""
    hypersurface(S::Spine2D)

Get the tropical curve which describes the spine.
"""
hypersurface(S::Spine2D) = S.curve

"""
    amoeba_approximation(S::Spine2D)

Get the approximated amoeba used to compute the spine.
"""
amoeba_approximation(S::Spine2D) = S.amoeba_approximation


function Base.show(io::IO, mime::MIME"text/html", S::Spine2D)
    if inIJulia() && plotsdefined()
        show(io, mime, Main.Plots.plot(S))
    else
        show(io, S)
    end
end
function Base.show(io::IO, S::Spine2D)
    nbounded = length(filter(isbounded, components_complement(S)))
    println(io, typeof(S), " with:")
    println(io, " * ", length(vertices(S.curve)), " vertices")
    println(io, " * ", length(halfrays(S.curve)), " half-rays")
    print(io, " * ", nbounded, " bounded cells")
end


@recipe function plot(S::Spine2D)
    xmin, xmax, ymin, ymax = limits(S.amoeba_approximation.grid)
    xlims --> (xmin, xmax)
    ylims --> (ymin, ymax)
    @series begin
        S.curve
    end
end
