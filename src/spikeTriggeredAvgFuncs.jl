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

function spikeTriggeredAvgNDimFilter(stimulus, response, windowSize, filter=[1], cond=1)
	spikeIndices = findall(x-> x == cond, sum(conv(response, filter), dims=2))
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
