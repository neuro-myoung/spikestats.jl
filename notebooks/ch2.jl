### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ cbadf790-3b3a-11eb-1057-317e8e2bab47
using Random, Distributions, Plots, DSP, Measures

# ╔═╡ b834f760-3e8a-11eb-07de-e9f27e984ca8
begin
	using MAT
	dat = matread("../data/c1p8.mat")
end;

# ╔═╡ 84bb19d0-3e8c-11eb-185c-c5153006e5af
md"# Exercise 1"

# ╔═╡ 0ac4a1e0-3b3b-11eb-07c9-051a8d640b9a
begin
	
	function whiteNoiseGenerator(T, interval, σ, seed=123)
		Random.seed!(seed)
		d = Normal(0, σ^2/interval)
		return rand(d, Int(T/interval))
	end
	
	function D(τ)
		return -cos(2π*(τ-0.02))/0.14 * exp(-τ/0.06)
	end
	
end

# ╔═╡ 28099d00-3b3b-11eb-0fa0-4d88f26e5869
begin
	s = whiteNoiseGenerator(10,0.01,sqrt(0.01));
	Ds = zeros(length(s));
	t = 0.01:0.01:10;
	for i in 1:length(t)
		Ds[i] = D(t[i]);
	end
	
	r0 = 50;
	rEst = 50 .+ conv(Ds,s)[1:1000];
	
corr = conv(rEst, reverse(s))[1:1000] ./ 10;
end

# ╔═╡ 547dc550-3b3b-11eb-0a74-e1e4fe59b5a3
begin
	plotly()
	
	l = @layout [a; b; c; d]
	
	p1 = plot(t, s, title="White Noise Stimulus",
		ylab="Stimulus Intensity (au)", xlab="Time (s)", legend=false)
	
	p2= scatter(t, Ds, legend=false, title="Filter")
	
	p3 = plot(t, rEst, title="Estimated Firing Rate",
		ylab="rEst (Hz)", xlab="Time (s)", legend=false)
	
	p4 = plot(t, corr, title="Estimated Firing Rate",
		ylab="rEst (Hz)", xlab="Time (s)", legend=false)
	
	plot(p1, p2, p3, p4, layout=l, size=(600,900), bottom_margin=10mm)
end

# ╔═╡ c8e44a00-3e8c-11eb-3599-85720a2de20a
md"# Exercise 2"

# ╔═╡ 1f3dd970-3e8d-11eb-3d6c-4f9625c97ee0
function spikeTriggeredAvg(stimulus, response, windowSize, filter=[1], cond=1)
	spikeIndices = findall(x-> x == cond, conv(response, filter))
	temp = zeros(Int(length(spikeIndices)), windowSize)
	for i in 1:1:length(spikeIndices)
		if spikeIndices[i][1] <= windowSize
			subArr = stimulus[1:spikeIndices[i][1]]
			prepend!(subArr, zeros(windowSize-spikeIndices[i][1]))
			temp[i,:] = subArr
		else
			temp[i,:] = stimulus[spikeIndices[i][1]-(windowSize-1):spikeIndices[i][1]]
		end
	end
	
	return mean(temp, dims=1), sqrt.(var(temp, dims=1)./size(temp)[1])
end

# ╔═╡ ffb81612-3e8c-11eb-0781-9142009ff055
begin
	m1, sem1 = spikeTriggeredAvg(dat["stim"], dat["rho"], 150);
	#rEst2 = 50 .+ conv(m1,dat["stim"])[1:1000]
end

# ╔═╡ 2d80aad0-3e8d-11eb-1f8c-e98354568725
plot((-0.298:0.002:0).*1000, m1[:], ribbon=sem1,fillalpha=.5, legend=false, 				xlabel="Time (ms)", ylabel="Stimulus Intensity (au)")

# ╔═╡ Cell order:
# ╟─84bb19d0-3e8c-11eb-185c-c5153006e5af
# ╠═cbadf790-3b3a-11eb-1057-317e8e2bab47
# ╠═0ac4a1e0-3b3b-11eb-07c9-051a8d640b9a
# ╠═28099d00-3b3b-11eb-0fa0-4d88f26e5869
# ╠═547dc550-3b3b-11eb-0a74-e1e4fe59b5a3
# ╟─c8e44a00-3e8c-11eb-3599-85720a2de20a
# ╠═b834f760-3e8a-11eb-07de-e9f27e984ca8
# ╠═1f3dd970-3e8d-11eb-3d6c-4f9625c97ee0
# ╠═ffb81612-3e8c-11eb-0781-9142009ff055
# ╠═2d80aad0-3e8d-11eb-1f8c-e98354568725
