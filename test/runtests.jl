using SpikeStats
using Test

testset = poissonSpikes(10,100,0.001)

println(names(SpikeStats))

@testset "SpikeStats.jl" begin
    @test typeof(testset) == Array{Float64,1}
    @test round(computeCV(diff(testset)), digits=2) == 0.93
    @test round(computeFano(testset,0.001), digits=2) == 0.90
end
