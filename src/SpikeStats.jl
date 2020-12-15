module SpikeStats

using Random
using StatsBase
using Statistics
using DataFrames
using DSP
using Distributions

include("poissonSpikeModels.jl")
include("computeCV.jl")
include("computeFano.jl")
include("acf.jl")
include("spikeTriggeredAvgFuncs.jl")
include("whiteNoiseGenerator.jl")
include("ccf.jl")

export poissonSpikes, computeCV, computeFano, poissonSpikesRefractory, poissonSpikesVarRate, acf, spikeTriggeredAvg, 
whiteNoiseGenerator, spikeTriggeredAvgNDimFilter, ccf, poissonSpikesStim
end
