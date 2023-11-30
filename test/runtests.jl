using ThomasAlgorithm
using Test

include("Poisson_test.jl")

@testset "ThomasAlgorithm.jl" begin
    @test Poisson_equation_test(0.1)
end