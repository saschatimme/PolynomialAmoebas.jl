@testset "Greedy" begin
    @polyvar x y z w
    f = x^2+y^2+1
    grid = Grid2D((-2, 2), 800)
    A = amoeba(f, grid=grid, alg=Greedy())
    @test any(A.data)

    I = imaginary_projection(f, grid=grid, alg=Greedy())
    @test any(I.data)

    f = x^2+y^2+z^2+1
    F = AmoebaFiber3D(f)
    grid = Grid3D(xlims=(-2, 2), ylims=(-2, 2), zlims=(-2, 2), res=(21, 21, 21))
    A3 = PolynomialAmoebas.amoeba(f, grid=grid, alg=Greedy())
    @test any(A3.data)

    I3 = imaginary_projection(f, grid=grid, alg=Greedy())
    @test any(I3.data)

    f = 1 + z * cis(0.8π) + z^3*w*cis(1.6π) + z*w^2*cis(0.4π) + w*cis(1.2π)
    C = coamoeba(f, alg=Greedy(), resolution=500)
    @test any(C.data)
end
