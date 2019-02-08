@testset "MembershipTest" begin
    @polyvar x y
    f = x^2 + y^2 + 1

    F = AmoebaFiber2D(f)

    result = membershiptest(F, (0.0, 0.0), PolynomialAmoebas.TorusStartValueGenerator{2}())
    @test result.successfull

    r2 = membershiptest(F, (0.1, 0.1), PolynomialAmoebas.TorusStartValueGenerator{2}(), (result.solution,), MembershipTestOptions())
    @test r2.successfull
    @test r2.tries == 1
end
