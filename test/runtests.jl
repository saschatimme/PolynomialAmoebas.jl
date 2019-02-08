using PolynomialAmoebas
using Test

@testset "PolynomialAmoebas" begin
    include("fibers.jl")
    include("membership_test.jl")
    include("newton_test.jl")

    include("newton_polygon_test.jl")
    include("tropical_test.jl")
    include("tropicalcurve_test.jl")

    include("simple_grid_test.jl")
    include("amoeba_archtrop_test.jl")
    include("contour_test.jl")
    include("greedy_grid_test.jl")
    include("spine_test.jl")
    include("polygonal_amoeba_test.jl")
end
