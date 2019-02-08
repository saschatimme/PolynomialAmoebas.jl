@testset "Newton polygon" begin

    @polyvar x y
    f = 4 + 4x + 4y + 4x*y + 2x^2 + 2y^2

    newtonpoly = newtonpolygon(f)
    @test lattices(newtonpoly) ==[[2, 0], [1, 1], [0, 2], [1, 0], [0, 1], [0, 0]]
    @test vertices(newtonpoly) == [[2, 0], [0, 2], [0, 0]]
    @test vertexindices(newtonpoly) == [1, 3, 6]
    @test facets(newtonpoly) == [[[2, 0], [0, 2]], [[0, 2], [0, 0]],[[0, 0], [2, 0]]]
    @test newtonpoly.subdivision == [[5, 2, 3], [4, 2, 1], [2, 5, 6, 4]]

    @test_throws AssertionError newtonpolygon(1 + x)

    @test_throws AssertionError vertices(newtonpolygon(x + y))

    @polyvar z w
    f = 1 + z * cis(0.8π) + z^3*w*cis(1.6π) + z*w^2*cis(0.4π) + w*cis(1.2π)
    n = newtonpolygon(archimedean_tropical_polynomial(f))
    @test length(n.subdivision) == 1
    @test length(n.subdivision[1]) == 5

    @polyvar x y
    f = x^2*y + y^2 + 3x^2*y^3 + y^4 + x^4*y^4
    newt = newtonpolygon(f)
    @test PolynomialAmoebas.on_boundary(newt, [2, 3]) == false
    @test PolynomialAmoebas.iscontained(newt, [2, 3])
end
