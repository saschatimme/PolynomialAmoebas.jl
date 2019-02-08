import StaticArrays: SVector

@testset "Newton" begin
    @polyvar x y z

    f = x^2*y + y^2 + 3x^2*y^3 + y^4 + x^4*y^4

    AF = AmoebaFiber2D(f, (0, 0))
    PolynomialAmoebas.update_fiber!(AF, (0, 1))
    @test PolynomialAmoebas.newton(AF, SVector(0.5π, 1.5π), 100, 1e-6).converged

    # underdetermined
    AF = AmoebaFiber3D(x^2+y^2+z^2+1, (0, 1, 0))
    PolynomialAmoebas.update_fiber!(AF, (0, 0, 0))
    v = SVector(3.97394, 2.84526, 1.69766)
    @test PolynomialAmoebas.newton(AF, v, 150, 1e-7).converged
end
