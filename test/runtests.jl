using SpikeStats
using Test

testset = poissonSpikes(10,100,0.001)

function rateFunc(t)
	return 100(1 + cos(2*π*t/0.025))
end

println(names(SpikeStats))

testset = poissonSpikes(10,100,0.001)

function rateFunc(t)
	return 100(1 + cos(2*π*t/0.025))
end

println(names(SpikeStats))

@testset "SpikeStats.jl" begin
    @test typeof(testset) == Array{Float64,1}
    @test round(mean(testset), digits=2) == 4.91
    @test typeof(poissonSpikesRefractory(20,100,0.02,0.001)[1]) == Array{Float64,2}
    @test round(mean(poissonSpikesRefractory(20,100,0.02,0.001)[1]),digits=2) == 10.02
    @test round(computeCV(diff(testset)), digits=2) == 1.07
    @test round(computeFano(testset,0.001), digits=2) == 0.90
    @test typeof(acf(testset, 0.001, 0.1)[2]) == Array{Float64,1}
    @test round(mean(acf(testset, 0.001, 0.1)[2]), digits=2) == 2.67
    @test round(mean(poissonSpikesVarRate(50,0.001,rateFunc)[1]), digits=2) == 25.11
    @test typeof(poissonSpikesVarRate(50,0.001,rateFunc)[1]) == Array{Float64,2}
end
