@testset "Fibers" begin
    @polyvar x y z
    f = x^2+y^2+1

    AF = AmoebaFiber2D(f)
    PolynomialAmoebas.update_fiber!(AF, (0.0, 0.0))

    @test_throws AssertionError AmoebaFiber2D(x^2+1)
    @test_throws AssertionError AmoebaFiber2D(x^2+z+y)

    CF = CoamoebaFiber2D(f)
    PolynomialAmoebas.update_fiber!(CF, (0.0, 0.0))
    @test_throws AssertionError CoamoebaFiber2D(x^2+1)
    @test_throws AssertionError CoamoebaFiber2D(x^2+z+y)

    CF2 = ContourFiber2D(f)
    PolynomialAmoebas.update_fiber!(CF2, 0.0)
    PolynomialAmoebas.fix_axis!(CF2, :w2)
    PolynomialAmoebas.update_fiber!(CF2, 1.0)
    @test_throws AssertionError ContourFiber2D(x^2+1)
    @test_throws AssertionError ContourFiber2D(x^2+z+y)

    I = ImaginaryFiber2D(f)
    @test_throws AssertionError ImaginaryFiber2D(x^2+1)

    g = x^2+y^2+z^2+1
    AG = AmoebaFiber3D(g)
    PolynomialAmoebas.update_fiber!(AG, (0.0, 0.0, 0.0))
    @test_throws AssertionError AmoebaFiber3D(x^2+1)
    @test_throws AssertionError AmoebaFiber3D(x^2+z)

    CG = CoamoebaFiber3D(g)
    PolynomialAmoebas.update_fiber!(CG, (0.0, 0.0, 1.0))
    @test_throws AssertionError CoamoebaFiber3D(x^2+1)
    @test_throws AssertionError CoamoebaFiber3D(x^2+z)

    IG = ImaginaryFiber3D(g)
    @test_throws AssertionError ImaginaryFiber3D(x^2+1)
end
