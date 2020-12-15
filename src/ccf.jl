function ccf(arr1, arr2, binWidth, T, demean=true)
	dtArr = (arr1 .- arr2')[:]
	histArr = fit(Histogram, filter(x->abs(x)<=T,dtArr), 
				  -T+0.5binWidth:binWidth:T-0.5binWidth)
	binLocs = -T+0.5binWidth:binWidth:T-0.5binWidth
	
	if demean == true
		weightsDemeaned = (histArr.weights 
							.- (sum(histArr.weights)*binWidth/(2T)))/(2T)
		return binLocs, weightsDemeaned
	else 
		return binLocs, histArr.weights
	end
end