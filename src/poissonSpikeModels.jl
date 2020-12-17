function poissonSpikes(T,freq,interval,seed=1)
	rng = MersenneTwister(seed)
	randArr = rand!(rng, zeros(Int(T/interval)))
	if isa(freq, Array)
		spikes = findall(x -> x == true, randArr .< freq .* interval) .* interval
	else
		spikes = findall(x -> x < freq * interval, randArr)/(1/interval)
	end
	
	return spikes
end

function poissonSpikesRefractory(T,freq0,tau,interval,seed=1)
	rng = MersenneTwister(seed)
	df = DataFrame(time = interval:interval:T, 
				   prob = rand!(rng, zeros(Int(T/interval))))
	df.freq = zeros(length(df.time))
	df.spikes = zeros(length(df.time))
	prevSpike = -10000
	for i in 1:size(df, 1)
		dt = df.time[i] - prevSpike
		freq = freq0 * (1 - exp(-1 * ( dt/tau)))
		df.freq[i] = freq
		if df.prob[i] < freq*interval
			df.spikes[i] = 1.0
			prevSpike = df.time[i]
		else
			df.spikes[i] = 0
		end
	end
		
	return df.time[df[:spikes] .== 1, :], df
end

function poissonSpikesVarRate(T, interval, func, seed=1)
	rng = MersenneTwister(seed)
	df = DataFrame(time = interval:interval:T, 
				   prob = rand!(rng, zeros(Int(T/interval))))
	df.freq = zeros(length(df.time))
	df.spikes = zeros(length(df.time))
	for i in 1:size(df, 1)
		freq = func(i*interval)
		df.freq[i] = freq
		if df.prob[i] < freq*interval
			df.spikes[i] = 1.0
			prevSpike = df.time[i]
		else
			df.spikes[i] = 0
		end
	end
		
	return df.time[df[:spikes] .== 1, :], df
end

function poissonSpikesStim(T, s, freq0,tau,interval,seed=1)
	rng = MersenneTwister(seed)
	df = DataFrame(time = interval:interval:T, 
					 prob = rand!(rng, zeros(Int(T/interval))))
	df.freq = zeros(length(df.time))
	df.spikes = zeros(length(df.time))
	prevSpike = 0
	for i in 1:size(df, 1)
		dt = df.time[i] - prevSpike
		df.freq[i] = abs.(s[i] .+ freq0) .- exp(-dt/tau)  
		if df.prob[i] < df.freq[i] * interval
			df.spikes[i] = 1.0
			prevSpike = df.time[i]
		else
			df.spikes[i] = 0
		end
	end

	return df.time[df[:spikes] .== 1, :], df
end