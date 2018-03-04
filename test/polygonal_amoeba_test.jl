@testset "PolygonalAmoeba" begin
    @polyvar x y
    f = x^2*y + y^2 + 3x^2*y^3 + y^4 + x^4*y^4

    A = amoeba(f, accuracy=0.01)
    @test 0 < accuracy(A) < 0.01
end
