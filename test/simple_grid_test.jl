@testset "SimpleGrid" begin
    @polyvar x y z w
    f = x^2+y^2+1
    F = AmoebaFiber2D(x^2+y^2+1)
    grid = Grid2D((-2, 2), 50)
    generator = PolynomialAmoebas.startvalue_generator(F)
    B = PolynomialAmoebas.simple_grid(F, grid, generator)
    @test any(B.data)
    B = PolynomialAmoebas.amoeba(f, grid=grid, alg=Simple())
    @test any(B.data)

    I = PolynomialAmoebas.imaginary_projection(f, grid=grid, alg=Simple())
    @test any(I.data)

    f = x^2+y^2+z^2+1
    F = AmoebaFiber3D(f)
    grid = Grid3D(xlims=(-2, 2), ylims=(-2, 2), zlims=(-2, 2), res=(21, 21, 21))
    generator = PolynomialAmoebas.startvalue_generator(F)
    B = PolynomialAmoebas.simple_grid(F, grid, generator)
    @test any(B.data)
    B = PolynomialAmoebas.amoeba(f, grid=grid, alg=Simple())
    @test any(B.data)

    f = 1 + z * cis(0.8π) + z^3*w*cis(1.6π) + z*w^2*cis(0.4π) + w*cis(1.2π)
    C = CoamoebaFiber2D(f)
    grid = Grid2D((0.0, 2π), 50)
    generator = PolynomialAmoebas.DomainStartValueGenerator2D((-5, 5), (-5, 5))
    B = simple_grid(C, grid, generator)
    @test any(B.data)
    B = coamoeba(f, alg=Simple())
    @test any(B.data)
end
