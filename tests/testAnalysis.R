library(evaluomeR)

data("rnaMetrics")
plotMetricsMinMax(rnaMetrics)
plotMetricsBoxplot(rnaMetrics)
plotMetricsCluster(rnaMetrics)
plotMetricsViolin(rnaMetrics)

stabilityData <- stabilityRange(data=rnaMetrics, k.range=c(3,4), bs=20, getImages = FALSE, seed=100)
qualityData <- qualityRange(data=rnaMetrics, k.range=c(3,4), getImages = FALSE, seed=100)

kOptTable <- getOptimalKValue(stabilityData, qualityData, k.range=c(3,4))
kOptTable


df = assay(rnaMetrics)
k.vector1=rep(5,length(colnames(df))-1)
k.vector2=rep(2,length(colnames(df))-1)
plotMetricsClusterComparison(df, k.vector1, k.vector2)
