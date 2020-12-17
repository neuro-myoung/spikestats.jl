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
include("spikeCorr.jl")
include("spikeTriggeredAvgFuncs.jl")
include("whiteNoiseGenerator.jl")

export poissonSpikes, computeCV, computeFano, poissonSpikesRefractory, poissonSpikesVarRate,spikeTriggeredAvg, spikeCorr, whiteNoiseGenerator, spikeTriggeredAvgNDimFilter, poissonSpikesStim
end
