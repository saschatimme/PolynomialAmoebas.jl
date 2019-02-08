@testset "Contour" begin
    @polyvar x y
    p1 = 1 + x + y + x * y + y^2 + x^2 * y

    C = PolynomialAmoebas.contour(p1, res=(100, 100))
    @test length(C.ps) != 0
end
