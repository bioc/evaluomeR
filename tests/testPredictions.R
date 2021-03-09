library(evaluomeR)

globalMetricDf = globalMetric(rnaMetrics, k.range = c(2,3), nrep=10, criterion="AIC", PCA=TRUE)
globalMetricDf
