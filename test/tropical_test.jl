@testset "Tropical" begin

    a = Tropical(3.0)
    @test Tropical(2.3) isa Tropical{Float64}
    @test Tropical(2) isa Tropical{Int64}
    @test inf(2.0) isa Tropical{Float64}


    @test isinf(inf(Int))

    b = Tropical(2.0)

    @test a + b == a
    @test a * b == Tropical(5.0)

    @test (Tropical(3.0) == Tropical(3.0)) == true
    @test (Tropical(3.1) == Tropical(3.0)) == false

    @test one(typeof(Tropical(2))) == Tropical(0)
    @test one(Tropical(3.0)) == Tropical(0.0)

    @test zero(typeof(Tropical(3))) == inf(Int)
    @test zero(Tropical(2.3)) == inf(Float64)

    @test inf(Tropical(3)) == inf(Int)


    @test Tropical(3) + Tropical(5.2) == Tropical(5.2)
    @test inf(Int) + Tropical(5.2) == Tropical(5.2)

    @test string(Tropical(5.2)) == "5.2"
    @test string(inf(Int)) == "-∞"
end
