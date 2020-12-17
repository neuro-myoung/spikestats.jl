
"""
		spikeCorr(arr1, binWidth, T)
		spikeCorr(arr1, arr2, binWidth, T)

Spike timing correlation function. If one array is given calculates autocorrelation if two arrays are given calculates cross-correlation.

# Arguments
arr1=1-dimensional array of spike timings\\
\\
binWidth=size of bin interval in seconds\\
\\
T=Total time overwhich to calulate correlation\\

# Output
corrArr=A 1-dimensional array of correlation coefficients

# Example
```julia-repl
julia> spikeCorr(arr1,0.002, 0.1)
Float64[0.089335, 0.01514, 0.0151317, 0.02228, ...]

```
```julia-repl
julia> spikeCorr(arr1, arr2, 0.002, 0.1)
Float64[0.0893333, 0.0207067, 0.00866, 0.00900833, ...]
```
"""
function spikeCorr(arr1, binWidth, T)
	corrArr = zeros(Int(T/binWidth))
	
	for i in 1:Int(T/binWidth)
		corrArr[i] = sum(fit(Histogram, arr1, (0:0.002:1200)).weights .*
						 fit(Histogram, arr1 .+ (i-1)*binWidth, 											(0:0.002:1200)).weights)/(1200/0.002)
	end
	return corrArr
end

function spikeCorr(arr1, arr2, binWidth, T)
	corrArr = zeros(Int(T/binWidth))
	
	for i in 1:Int(T/binWidth)
		corrArr[i] = sum(fit(Histogram, arr1, (0:0.002:1200)).weights .*
						 fit(Histogram, arr2 .+ (i-1)*binWidth, 											(0:0.002:1200)).weights)/(1200/0.002)
	end
	return corrArr
end