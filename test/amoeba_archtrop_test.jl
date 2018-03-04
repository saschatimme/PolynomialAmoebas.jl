@testset "ArchTropAmoeba" begin
    @polyvar x y
    f = x^2+y^2+1

    grid = Grid2D((-2, 2), 30)
    B = amoeba(f, grid=grid, alg=ArchTrop())
    @test any(B.data)
end
