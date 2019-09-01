@testset "Spine2D" begin
    @polyvar x y
    f = x^2*y + y^2 + 3x^2*y^3 + y^4 + x^4*y^4

    S = spine(f)
    @test length(components_complement(S)) == 5
    @test sprint(show, S) ==
        """Spine2D{DynamicPolynomials.Polynomial{true,Tropical{Float64}}} with:
            $(" ")* 4 vertices
            $(" ")* 4 half-rays
            $(" ")* 1 bounded cells"""
    S = spine(f, minimal_component_size=0.05)
    @test length(components_complement(S)) == 5

    S = spine(f, domain=(-5, 5, -5, 5), minimal_component_size=0.01)
    @test length(components_complement(S)) == 5
end
