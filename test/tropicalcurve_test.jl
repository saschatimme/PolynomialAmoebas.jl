@testset "Tropical Curve" begin
    @polyvar x y

    g = Tropical(0) + x^3 + y^3 + log(2)*x*y
    t = TropicalCurve(g)
    n3 = [0.0, -log(2)]
    n2 = [-log(2), 0.0]
    n1 = [log(2), log(2)]
    @test vertices(t)[1] ≈ n1
    @test vertices(t)[2] ≈ n3
    @test vertices(t)[3] ≈ n2

    @test segments(t)[1][1] ≈ n1
    @test segments(t)[1][2] ≈ n3
    @test segments(t)[2][1] ≈ n1
    @test segments(t)[2][2] ≈ n2
    @test segments(t)[3][1] ≈ n3
    @test segments(t)[3][2] ≈ n2

    @test length(halfrays(t)) == 3


    # c = ArchTrop(1 + x + y + x * y + y^2 + x^2 * y).curve
    # @test length(halfrays(c)) == 4
    # @test length(vertices(c)) == 1
end
