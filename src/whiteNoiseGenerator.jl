function whiteNoiseGenerator(T, interval, sigma, seed=123)
	Random.seed!(seed)
	d = Normal(0,sigma^2/interval)
	return rand(d, Int(T/interval))
end