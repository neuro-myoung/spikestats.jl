function computeFano(arr, interval)
  arrHist = fit(Histogram, arr, 0:interval:Int(round(maximum(arr))))
  ff = (mean(arrHist.weights.^2)-mean(arrHist.weights)^2)/mean(arrHist.weights)
  return ff
end