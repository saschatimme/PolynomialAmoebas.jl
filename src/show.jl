
Base.show(io::IO, F::AbstractFiber) = print(io, typeof(F), "()")

function Base.show(io::IO, mime::MIME"text/html", NP::NewtonPolygon)
    if inIJulia() && plotsdefined()
        show(io, mime, Main.Plots.plot(NP))
    else
        show(io, NP)
    end
end
function Base.show(io::IO, NP::NewtonPolygon)
    #println(typeof(A), ":")
    println(io, "Newton polygon")
    println(io, "* lattices: ", lattices(NP))
    println(io, "* vertices: ", vertices(NP))
end

function Base.show(io::IO, B::Bitmap2D)
    n, m = size(B)
    print(io, "Bitmap2D of size ", n, "×", m)
end
function Base.show(io::IO, B::Bitmap3D)
    n1, n2, n3 = size(B)
    print(io, "Bitmap2D of size ", n1, "×", n2, "×", n3)
end
function Base.show(io::IO, mime::MIME"text/html", M::AbstractBitmap)
    if inIJulia() && plotsdefined()
        show(io, mime, Main.Plots.plot(M))
    else
        show(io, M)
    end
end
