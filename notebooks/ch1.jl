### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ d6111790-0cf4-11eb-1187-a9b6c6018291
using SpikeStats, Plots, StatsPlots, Measures, Random

# ╔═╡ 74a7f0b0-1257-11eb-29cd-2dcc523aaeea
begin
	using StatsBase
	using DSP
	
	l6 = @layout [a;b;c]
	
	stim = whiteNoiseGenerator(100,0.001,0.005)
	p6a = plot(0.001:0.001:100, stim, linetype=:steppre, legend=false, 
		xlabel="Time (s)", title="Signal")
	
	x5 = autocor(stim; demean=true)
	p6b = bar(x5, color=:black, fill=:black, bar_width=1, legend=false,
		title="Autocorrelation", xlabel="Time (ms)")
	
	x6 = periodogram(stim, fs=1000, window=nothing)
	p6c = plot(x6.freq,x6.power, xlabel="Frequency (Hz)", ylabel="Power",
		legend=false, title="Power Spectral Density")
	
	plot(p6a, p6b, p6c, layout=l6, size=(600,900), bottom_margin=15mm)
end

# ╔═╡ 60f14d52-125c-11eb-3d0a-719cb60752dd
begin
	using MAT
	dat = matread("../data/c1p8.mat")
end;

# ╔═╡ 14e29010-0ef9-11eb-06c8-f50d6ec0c59d
md"# Exercise 1"

# ╔═╡ 2237f990-0f0c-11eb-3983-1bf95b4a958f
md"## Poisson Spike Raster"

# ╔═╡ a77e632e-0ef0-11eb-0b65-5543bc1d2e94
@bind lastSpike html"<input type=range min=10 max=1000 step=10>"

# ╔═╡ b7244700-0ef0-11eb-30e5-31fd6ab4e2af
md"Showing the first $lastSpike spikes."

# ╔═╡ 0d5f3500-0cfb-11eb-2525-3d5a35a84ce3
begin
	model1 = poissonSpikes(20,100,0.001);
	
	plotly()
	scatter(model1[1:lastSpike], fill(1, lastSpike), marker=:vline, ylim=[0,2], 				yaxis=false, grid=false, legend=false, xlab="Time (s)")
end

# ╔═╡ 2dce4800-3a8f-11eb-2bdd-cd663793c4cd
begin
	v = [4.11e-21, 8e-21, 1.6e-19, 1e-15, 100e-15]
	plotly()
	scatter(v, fill(1, 5), marker=:vline, ylim=[0,2], yaxis=false, grid=false, legend=false, xlab="", xaxis=:log,
	xrotation=-90, xtickfont=font(12))
end

# ╔═╡ 93c90050-3a8f-11eb-15c4-b9f6fe0e2365
v

# ╔═╡ 986da92e-3a8f-11eb-1f4c-6517c23699dc


# ╔═╡ dbe518b0-0cf4-11eb-17d7-c343bd190040
isis = diff(model1);

# ╔═╡ d03cd250-0ef7-11eb-243b-2586ebc58a54
md"Slider for interval from 1 to 100 ms"

# ╔═╡ f9f07080-0ef6-11eb-30fd-254c1c5b6e2d
@bind interval html"<input type=range min=0.001 max=0.1 step=0.001>"

# ╔═╡ dc49f472-0ef6-11eb-24b0-5b84150e719f
begin
	cv = computeCV(isis)
	md"Coefficient of Variation: $cv"
end

# ╔═╡ 1fbd62ee-2ab3-11eb-28ee-8b0a384b3cfa
begin
	ff = computeFano(model1[:], interval)
	md"FanoFactor: $ff"
end

# ╔═╡ 347acf60-0f0c-11eb-37a3-69b6207a4b61
md"## Poisson Interspike Interval Histogram"

# ╔═╡ 68741b50-0ef8-11eb-16a4-d71c74789eb0
histogram(isis, bins = :scott, legend=false, xlab="Time(s)", ylab="Counts")

# ╔═╡ 607332a0-0ef9-11eb-3377-615569385b52
md"# Exercise 2"

# ╔═╡ d80e1d80-0f0c-11eb-2250-df712fba9635
md"## Ex Simulation with tau=20 ms"

# ╔═╡ 6dab8920-0f0a-11eb-3dae-8f35e4db9672
md"View first n seconds of simulation"

# ╔═╡ 529506c0-0f0a-11eb-22d6-d7852b19a4aa
@bind lastSpike2 html"<input type=range min=0 max=100 step=0.1>"

# ╔═╡ 4ced0232-0f0b-11eb-1265-7df1435f19fd
md"Showing first $lastSpike2 seconds."

# ╔═╡ 46bf04c0-0f0c-11eb-0e38-6d4944aaa71f
md"## Poisson with Refractory Period"

# ╔═╡ 5f5038a0-0f08-11eb-1518-ddd6e1914a16
begin
	model2, model2df = poissonSpikesRefractory(20,100,0.02,0.001);
	
	l2a = @layout [a ; b]
	subDf = model2df[model2df.time .<= lastSpike2, :]
	subSpikeTimes = model2[model2 .<= lastSpike2]


	p1a = @df subDf plot(:time, :freq, grid=false, legend=false, 
		ylab="Effective Frequency (Hz)")
	p1b = 	scatter(subSpikeTimes, fill(1, size(subSpikeTimes,1)), marker=:vline, 				ylim=[0,2], yaxis=false, grid=false, legend=false, xlab="Time (s)")
	
	plot(p1a, p1b, layout = l2a)

	
	
end

# ╔═╡ 4772a080-0f10-11eb-3032-91f1ef745f01
begin
	taus =0.001:0.001:0.02
	spikesArr = fill(Float64[], 1, length(taus))
	isisArr =  fill(Float64[], 1, length(taus))
	cvArr = zeros(length(taus))

	for i = 1:length(taus)
		model2Refract,model2RefractDf = poissonSpikesRefractory(20,100,taus[i],0.001);
		spikesArr[i] = reduce(vcat, 
			model2RefractDf.time[model2RefractDf[:spikes] .== 1, :])
		isisArr[i] =  diff(reduce(vcat,spikesArr[i]))
		cvArr[i] = computeCV(isisArr[i])
	end
	
	l = @layout [a ; b c d]
	
	p1 = plot(taus, cvArr, ylabel="CV", xlabel="tau (s)",bottom_margin = 10mm)
	
	p2 = histogram(isisArr[1], bins = :scott, legend=false, ylab="Counts", 							   title="tau=1 ms")

	p3 = histogram(isisArr[10], bins = :scott, legend=false, xlab="Time(s)", 						   ylab="Counts", title="tau=10 ms")

	p4 = histogram(isisArr[20], bins = :scott, legend=false, ylab="Counts", 	 	 	                title="tau=20 ms")
	
	plot(p1, p2, p3, p4, layout = l)
end

# ╔═╡ 19473900-0f15-11eb-01b1-eb4e5b7aa5c6
@bind interval2 html"<input type=range min=0.001 max=0.1 step=0.001>"

# ╔═╡ 86b857d0-0f15-11eb-33d7-0bffb894f7a9
begin
	model2FF, model2FFDf = poissonSpikesRefractory(100,100,0.01,0.001);
	ff2 = computeFano(model2FF[:], interval2)
	md"Fano Factor with an interval of $interval2 seconds: $ff2"
end

# ╔═╡ f1c55270-0f16-11eb-1fae-5fd10c13576c
md"# Exercise 3"

# ╔═╡ a76b7030-124d-11eb-3e16-29c6179d954a
md"### Constant Firing Rate Poisson Process"

# ╔═╡ 6b977a90-0f28-11eb-30f1-9788f179ebf2
begin
	x,y = acf(model1, 0.001, 0.1)
	x2,y2 = acf(model2, 0.001, 0.1)
	
	function rateFunc(t)
		return 100(1 + cos(2*π*t/0.025))
	end
	
	model3, model3df = poissonSpikesVarRate(50,0.001,rateFunc);
	
	x3,y3 = acf(model3, 0.001, 0.1)
	
	
	l3 = @layout [a; b; c]
	p3a = bar(x.*1000, y,color=:black, fill=:black, bar_width=2, legend=false, 
			 ylabel="Counts (mean subtracted)", size=(600,300), title = "Poisson")
	p3b = bar(x2.*1000, y2,color=:black, fill=:black, bar_width=2, legend=false, 					ylabel="Counts (mean subtracted)", size=(600,300), 
				title = "Poisson + Refractory")
	p3c = bar(x3.*1000, y3,color=:black, fill=:black, bar_width=2, legend=false, 					xlabel="Time (ms)", ylabel="Counts (mean subtracted)", 
			title="Poisson Time-variant",size=(600,300))
	
	plot(p3a, p3b, p3c, layout = l3, size=(600,900), bottom_margin=10mm)
end

# ╔═╡ 58bef230-124e-11eb-094c-6158ea264c00
md"# Exercise 4"

# ╔═╡ 5434f770-1a08-11eb-0df4-4bcd9d419226
@bind window html"<input type=range min=0.05 max=2 step=0.05>"

# ╔═╡ 973144c0-1a26-11eb-3a88-c325b49710d9
md"Window: $window s"

# ╔═╡ 97c2fd70-1a26-11eb-2b39-b108b73024ed
@bind tauGuess html"<input type=range min=0.001 max=0.1 step=0.001>"

# ╔═╡ 4e93ca80-1a26-11eb-19df-99f51df8c143
begin
	
	function rateFunc2(t)
		return 100(1 + cos(2*π*t/0.3))
	end
		
	spikes5, spikes5df = poissonSpikesVarRate(10,0.001,rateFunc2);
	subDf5 = spikes5df[spikes5df.time .<= window, :];
	
	using DataFrames
	
	function poissonSpikesCustom(T,freq0,tau,interval,seed=1)
		rng = MersenneTwister(seed)
		df = DataFrame(time = interval:interval:T, 
				   	prob = rand!(rng, zeros(Int(T/interval))))
		df.freq = zeros(length(df.time))
		df.spikes = zeros(length(df.time))
		prevSpike = 0
		fSet= freq0
		for i in 1:size(df, 1)
			dt = df.time[i] - prevSpike
			df.freq[i] = fSet * exp(-dt/tau)  
			if df.prob[i] < fSet * interval
				df.spikes[i] = 1.0
				fSet = df.freq[i] + 1/tau
				prevSpike = df.time[i]
			else
				df.spikes[i] = 0
			end
		end
		
		return df.time[df[:spikes] .== 1, :], df
	end
	
	function mse(x, y)
		return (1/(length(x))) * sum((y-x).^2)
	end

	spikes6, spikes6df = poissonSpikesCustom(50, 100, tauGuess, 0.001)
	subDf6 = spikes6df[spikes6df.time .<= window, :];
	
	l4 = @layout [a b; c d]
	
	p4a = @df subDf5 plot(:time, :freq, grid=false, legend=false, 
				    ylab="Effective Frequency (Hz)", title="Ground Truth")
	p4b = scatter(spikes5[spikes5 .<= window], 
		fill(1, length(spikes5[spikes5 .<= window])), marker=:vline, 
		ylim=[0,2], yaxis=false, grid=false, legend=false, xlab="Time (s)")
	p4c = @df subDf6 plot(:time, :freq, grid=false, legend=false, 
				    ylab="Effective Frequency (Hz)", 
					title="Best Fit Estimate: $(tauGuess*1000) ms")
	p4d = scatter(spikes6[spikes6 .<= window], fill(1, 
			length(spikes6[spikes6 .<= window])), marker=:vline, 
		ylim=[0,2], yaxis=false, grid=false, legend=false, xlab="Time (s)")

	plot(p4a, p4c, p4b, p4d, layout = l4)
	
	
end

# ╔═╡ bb425ca0-1a26-11eb-353e-9b38f36dd536
md"tau: $(tauGuess*1000) ms"

# ╔═╡ 4936bd20-1fc3-11eb-030c-2deb51fb27a0
	begin
	
		z = zeros(100)
		for i in 1:1:100
			model4b, spikes6dfopt = poissonSpikesCustom(10, 100, i/1000, 0.001)
			z[i] = mse(spikes6dfopt.freq, spikes5df.freq)
		end
		
		tauVal, tauValInd = findmin(z)
	
		spikes6opt, spikes6dfopt = poissonSpikesCustom(10, 100, tauValInd/1000, 0.001)
		subDf6opt = spikes6dfopt[spikes6dfopt.time .<= window, :];
		
		l4b = @layout [a b; c d]
		
		p4e = @df subDf5 plot(:time, :freq, grid=false, legend=false, 
							ylab="Effective Frequency (Hz)", title="Ground Truth")
		p4f = scatter(spikes5[spikes5 .<= window], 
			fill(1, length(spikes5[spikes5 .<= window])), marker=:vline, 
			ylim=[0,2], yaxis=false, grid=false, legend=false, xlab="Time (s)")
		p4g = @df subDf6opt plot(:time, :freq, grid=false, legend=false, 
							ylab="Effective Frequency (Hz)", 
						title="Best Fit Estimate: $tauValInd ms")
		p4h = scatter(spikes6opt[spikes6opt .<= window], fill(1, 
				length(spikes6opt[spikes6opt .<= window])), marker=:vline, 
			ylim=[0,2], yaxis=false, grid=false, legend=false, xlab="Time (s)")
	
		plot(p4e, p4g, p4f, p4h, layout = l4b)
		
		
	end

# ╔═╡ 4b1229d0-124f-11eb-2596-375c8fff3e21
md"# Exercise 5"

# ╔═╡ 6a221230-1255-11eb-3aad-9511c2120b89
md"# Exercise 6"

# ╔═╡ 0d15e4a0-1259-11eb-0cf8-3f717eac7919
md"## White-noise simulation"

# ╔═╡ 915eb8c0-125b-11eb-3a78-4f41ed115fc0
md"# Exercise 7"

# ╔═╡ 894c56a0-2a94-11eb-3905-45bc52b71a27
WNStim = 1000*whiteNoiseGenerator(100,0.001,0.005)


# ╔═╡ a2b6bb60-2aaf-11eb-2a8c-599a9b55cbcd
plot(1:length(WNStim), abs.(WNStim .+ 90))

# ╔═╡ 975d9a70-125b-11eb-2204-bf4954f75ddd
md"# Exercise 8"

# ╔═╡ 2c2a2240-12ed-11eb-33ec-7b267c12df04
begin
	m1, sem1 = spikeTriggeredAvg(dat["stim"], dat["rho"], 150)
	plot((-0.298:0.002:0).*1000, m1[:], ribbon=sem1,fillalpha=.5, legend=false, 				xlabel="Time (ms)", ylabel="Stimulus Intensity (au)")
end

# ╔═╡ ed6f0d10-2ab4-11eb-0821-8fba6628a523
begin
	
	test,testdf = poissonSpikesStim(100, WNStim, 0.09, 0.02, 0.001);

	x7,y7 = ccf(WNStim, test, 0.001, 0.1)
	m7, sem7 = spikeTriggeredAvg(WNStim, testdf[:spikes], 150)
	
	l7 = @layout [a ; b]
	
	p7a = bar(x7.*1000, y7,color=:black, fill=:black, bar_width=2, legend=false, 
			 ylabel="Counts", size=(600,300), title = "Poisson")
	p7b = plot((-0.298:0.002:0).*1000, m7[:], ribbon=sem1,fillalpha=.5, legend=false, 				xlabel="Time (ms)", ylabel="Stimulus Intensity (au)")
	
	plot(p7a, p7b, layout=l7, size=(600,600))
end

# ╔═╡ 9d8b2d40-125b-11eb-2e1d-85d87aa9427d
md"# Exercise 9"

# ╔═╡ bbbac8f0-127b-11eb-022a-89bc3faef815
md"## Paired Spike Triggered Average"

# ╔═╡ 7c96ff40-127b-11eb-0fa2-93a150c3439b
@bind spikeDelay html"<input type=range min=2 max=100 step=2>"

# ╔═╡ d9589270-127b-11eb-3490-af073231ea5a
md"$spikeDelay millisecond delay"

# ╔═╡ b951be90-1283-11eb-15b0-6bfca37b3183
begin
	#single spikes
	m2, sem2 = spikeTriggeredAvg(dat["stim"], dat["rho"], 150)
	
	#pairs w/ delay
	gap = Int(spikeDelay/2)
	spikeFilter = vcat([1],zeros(gap-1),[1])
	mPair, semPair = spikeTriggeredAvg(dat["stim"], dat["rho"],  150, spikeFilter, 2)
	
	#single spikes delayed
	delayFilter = vcat(zeros(gap),[1]) 
	mDelay, semDelay = spikeTriggeredAvg(dat["stim"], dat["rho"], 150, delayFilter, 1)
	
	l2 = @layout [a ; b ; c]
	
	s1 = plot((-0.298:0.002:0).*1000, m2[:] + mDelay[1:150], ribbon=sem2,fillalpha=.5, 				legend=false, xlabel="Time (ms)", ylabel="Stimulus Intensity (au)",
			 size=(100,200), bottom_margin = 10mm)
	s2 = plot((-0.298:0.002:0).*1000, mPair[:], ribbon=semPair,fillalpha=.5, 					legend=false, xlabel="Time (ms)", ylabel="Stimulus Intensity (au)",
			size=(100,200), bottom_margin = 10mm)
	s3 = plot((-0.298:0.002:0).*1000, mPair[:] - (m1[:] + mDelay[1:150]), 						legend=false, xlabel="Time (ms)", ylabel="Stimulus Intensity (au)",
			size=(100,200))
	
	plot(s1, s2, s3, layout=l2, size=(600,600))
	
end

# ╔═╡ 22ba00f0-2a8e-11eb-3dcf-9f67664df0fb
begin
	meanSqDif = zeros(50)
	
	for i in 1:50
	
		gapo = Int((i*2)/2)
		spikeFiltero = vcat([1],zeros(gapo-1),[1])
		mPairo, semPairo = spikeTriggeredAvg(dat["stim"], dat["rho"],  150, 
			spikeFiltero, 2)
	
		delayFiltero = vcat(zeros(gapo),[1]) 
		mDelayo, semDelayo = spikeTriggeredAvg(dat["stim"], dat["rho"], 150,
			delayFiltero, 1)
	
		meanSqDif[i] = mse(mPairo, m2+mDelayo)
	
	end

	
	s4 = plot(2:2:100, meanSqDif, ylabel="Mean Sq. Dif.", xlabel="Gap Size (ms)", 
		size=(600,600), legend=false)
end

# ╔═╡ e6a2fab0-2a90-11eb-0e31-590e77e30ddb
md"# Exercise 10"

# ╔═╡ d36a61f0-2a8a-11eb-0bcc-43f097b53427
md"## Spike Delay with or without gap"

# ╔═╡ 62c95130-2a90-11eb-1c66-47619b9b794d
@bind spikeDelay2 html"<input type=range min=2 max=100 step=2>"

# ╔═╡ 74227960-2a87-11eb-2196-61cab3cecbc8
begin
	
	#pairs w/ delay
	gap2 = Int(spikeDelay2/2)
	spikeFilter3 = vcat([1],zeros(gap2-1),[1])
	mPair3, semPair3 = spikeTriggeredAvg(dat["stim"], dat["rho"],  150, spikeFilter3,
		2)
	
	spikeFilter2 = 2 .* vcat([1],ones(gap2-1),[1])
	mPair2, semPair2 = spikeTriggeredAvgNDimFilter(dat["stim"], dat["rho"], 
		150, hcat(spikeFilter3, spikeFilter2), 6)
	
	s5 = plot((-0.298:0.002:0).*1000, mPair2[:], ribbon=semPair2,fillalpha=.5,
		xlabel="Time (ms)", ylabel="Stimulus Intensity (au)",
			size=(600,400), bottom_margin = 10mm, legend=:topleft,
	label="Without spike gap")
	
	plot!((-0.298:0.002:0).*1000, mPair3[:], color="red", label="With spike gap")
end

# ╔═╡ Cell order:
# ╠═d6111790-0cf4-11eb-1187-a9b6c6018291
# ╟─14e29010-0ef9-11eb-06c8-f50d6ec0c59d
# ╟─2237f990-0f0c-11eb-3983-1bf95b4a958f
# ╟─a77e632e-0ef0-11eb-0b65-5543bc1d2e94
# ╟─b7244700-0ef0-11eb-30e5-31fd6ab4e2af
# ╠═0d5f3500-0cfb-11eb-2525-3d5a35a84ce3
# ╠═2dce4800-3a8f-11eb-2bdd-cd663793c4cd
# ╠═93c90050-3a8f-11eb-15c4-b9f6fe0e2365
# ╠═986da92e-3a8f-11eb-1f4c-6517c23699dc
# ╠═dbe518b0-0cf4-11eb-17d7-c343bd190040
# ╟─d03cd250-0ef7-11eb-243b-2586ebc58a54
# ╟─f9f07080-0ef6-11eb-30fd-254c1c5b6e2d
# ╟─dc49f472-0ef6-11eb-24b0-5b84150e719f
# ╟─1fbd62ee-2ab3-11eb-28ee-8b0a384b3cfa
# ╟─347acf60-0f0c-11eb-37a3-69b6207a4b61
# ╠═68741b50-0ef8-11eb-16a4-d71c74789eb0
# ╟─607332a0-0ef9-11eb-3377-615569385b52
# ╟─d80e1d80-0f0c-11eb-2250-df712fba9635
# ╟─6dab8920-0f0a-11eb-3dae-8f35e4db9672
# ╟─529506c0-0f0a-11eb-22d6-d7852b19a4aa
# ╟─4ced0232-0f0b-11eb-1265-7df1435f19fd
# ╟─46bf04c0-0f0c-11eb-0e38-6d4944aaa71f
# ╠═5f5038a0-0f08-11eb-1518-ddd6e1914a16
# ╠═4772a080-0f10-11eb-3032-91f1ef745f01
# ╟─19473900-0f15-11eb-01b1-eb4e5b7aa5c6
# ╠═86b857d0-0f15-11eb-33d7-0bffb894f7a9
# ╟─f1c55270-0f16-11eb-1fae-5fd10c13576c
# ╟─a76b7030-124d-11eb-3e16-29c6179d954a
# ╠═6b977a90-0f28-11eb-30f1-9788f179ebf2
# ╟─58bef230-124e-11eb-094c-6158ea264c00
# ╟─5434f770-1a08-11eb-0df4-4bcd9d419226
# ╟─973144c0-1a26-11eb-3a88-c325b49710d9
# ╟─97c2fd70-1a26-11eb-2b39-b108b73024ed
# ╟─bb425ca0-1a26-11eb-353e-9b38f36dd536
# ╠═4e93ca80-1a26-11eb-19df-99f51df8c143
# ╠═4936bd20-1fc3-11eb-030c-2deb51fb27a0
# ╟─4b1229d0-124f-11eb-2596-375c8fff3e21
# ╟─6a221230-1255-11eb-3aad-9511c2120b89
# ╟─0d15e4a0-1259-11eb-0cf8-3f717eac7919
# ╠═74a7f0b0-1257-11eb-29cd-2dcc523aaeea
# ╟─915eb8c0-125b-11eb-3a78-4f41ed115fc0
# ╠═894c56a0-2a94-11eb-3905-45bc52b71a27
# ╠═a2b6bb60-2aaf-11eb-2a8c-599a9b55cbcd
# ╠═ed6f0d10-2ab4-11eb-0821-8fba6628a523
# ╟─975d9a70-125b-11eb-2204-bf4954f75ddd
# ╠═60f14d52-125c-11eb-3d0a-719cb60752dd
# ╠═2c2a2240-12ed-11eb-33ec-7b267c12df04
# ╟─9d8b2d40-125b-11eb-2e1d-85d87aa9427d
# ╟─bbbac8f0-127b-11eb-022a-89bc3faef815
# ╟─7c96ff40-127b-11eb-0fa2-93a150c3439b
# ╟─d9589270-127b-11eb-3490-af073231ea5a
# ╠═b951be90-1283-11eb-15b0-6bfca37b3183
# ╠═22ba00f0-2a8e-11eb-3dcf-9f67664df0fb
# ╟─e6a2fab0-2a90-11eb-0e31-590e77e30ddb
# ╟─d36a61f0-2a8a-11eb-0bcc-43f097b53427
# ╟─62c95130-2a90-11eb-1c66-47619b9b794d
# ╠═74227960-2a87-11eb-2196-61cab3cecbc8
