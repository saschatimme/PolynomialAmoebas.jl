module PolynomialAmoebas

    using StaticArrays
    using RecipesBase
    using LinearAlgebra

    import MultivariatePolynomials
    const MP = MultivariatePolynomials
    import StaticPolynomials
    const SP = StaticPolynomials
    import CDDLib, Polyhedra

    import DynamicPolynomials: @polyvar

    export @polyvar

    import Colors
    import PlotUtils
    import Contour
    import DataStructures: PriorityQueue, enqueue!, dequeue!, peek, dequeue_pair!
    import Statistics: mean

    include("types.jl")
    include("utilities.jl")
    include("show.jl")

    include("fibers/fiber_helpers.jl")
    include("fibers/amoeba_fiber_2d.jl")
    include("fibers/amoeba_fiber_3d.jl")
    include("fibers/coamoeba_fiber_2d.jl")
    include("fibers/coamoeba_fiber_3d.jl")
    include("fibers/contour_fiber_2d.jl")
    include("fibers/imaginary_fiber_2d.jl")
    include("fibers/imaginary_fiber_3d.jl")

    include("newton.jl")

    include("membership_test_startvalues.jl")
    include("membership_test.jl")

    include("newtonpolygon.jl")

    include("grid_2d.jl")
    include("grid_3d.jl")
    include("bitmap_2d.jl")
    include("bitmap_3d.jl")
    include("bitmap.jl")

    include("tropical.jl")
    include("tropicalcurve.jl")
    include("archtrop.jl")

    include("simple_grid.jl")
    include("greedy_grid.jl")

    include("contour_2d.jl")
    include("contour.jl")

    include("amoeba.jl")
    include("coamoeba.jl")
    include("imaginary.jl")


    include("spine/components_complement.jl")
    include("spine/ronkin.jl")
    include("spine.jl")

    include("polygonal_amoeba/covering.jl")
    include("polygonal_amoeba/covering_fitter.jl")
    include("polygonal_amoeba.jl")
end
