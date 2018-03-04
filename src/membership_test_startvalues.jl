abstract type StartValueGenerator end


const two_pi = convert(Float64, 2π)

struct TorusStartValueGenerator{N} <: StartValueGenerator
end

startvalue(::TorusStartValueGenerator{2}, w=nothing) = SVector(rand() * two_pi, rand() * two_pi)
startvalue(::TorusStartValueGenerator{3}, w=nothing) = SVector(rand() * two_pi, rand() * two_pi, rand() * two_pi)


struct DomainStartValueGenerator2D <: StartValueGenerator
    xmin::Float64
    ymin::Float64

    dx::Float64
    dy::Float64
end

function DomainStartValueGenerator2D(lims::Tuple{<:Real, <:Real, <:Real, <:Real})
    xmin, xmax, ymin, ymax = convert.(Float64, lims)
    DomainStartValueGenerator2D(xmin, ymin, xmax - xmin, ymax - ymin)
end
DomainStartValueGenerator2D(xlims, ylims) = DomainStartValueGenerator2D((xlims..., ylims...))

# function DomainStartValueGenerator2D(f::MP.AbstractPolynomial; apply_exp=false, kwargs...)
#     xmin, xmax, ymin, ymax = amoeba_carcase_domain_heuristic(f; kwargs...)
#     DomainStartValueGenerator2D(xmin, ymin, xmax - xmin, xmax - ymin)
# end

function startvalue(c::DomainStartValueGenerator2D, y=nothing)
    w1 = rand() * c.dx + c.xmin
    w2 = rand() * c.dy + c.ymin

    SVector(w1, w2)
end


struct DomainStartValueGenerator3D <: StartValueGenerator
    xmin::Float64
    ymin::Float64
    zmin::Float64

    dx::Float64
    dy::Float64
    dz::Float64
end

function DomainStartValueGenerator3D(lims::Tuple{<:Real, <:Real, <:Real, <:Real, <:Real, <:Real})
    xmin, xmax, ymin, ymax, zmin, zmax = convert.(Float64, lims)
    DomainStartValueGenerator3D(
        xmin, ymin, zmin,
        xmax - xmin, ymax - ymin, zmax - zmin)
end
DomainStartValueGenerator3D(xlims, ylims, zlims) = DomainStartValueGenerator3D(xlims..., ylims..., zlims...)

function startvalue(c::DomainStartValueGenerator3D, w=nothing)
    w1 = rand() * c.dx + c.xmin
    w2 = rand() * c.dy + c.ymin
    w3 = rand() * c.dz + c.zmin

    SVector(w1, w2, w3)
end

startvalue_generator(::AmoebaFiber2D) = TorusStartValueGenerator{2}()
startvalue_generator(::AmoebaFiber3D) = TorusStartValueGenerator{3}()
