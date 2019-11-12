library(evaluomeR)

data("rnaMetrics")
plotMetricsMinMax(rnaMetrics)
plotMetricsBoxplot(rnaMetrics)
plotMetricsCluster(rnaMetrics)
plotMetricsViolin(rnaMetrics)

stabilityData <- stabilityRange(data=rnaMetrics, k.range=c(2,4),
                                bs=20, getImages = FALSE, seed = 100)
qualityData <- qualityRange(data=rnaMetrics, k.range=c(2,4),
                            getImages = FALSE, seed = 100)

kOptTable <- getOptimalKValue(stabilityData, qualityData, k.range=c(3,4))
kOptTable
