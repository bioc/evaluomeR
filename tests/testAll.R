library(evaluomeR)

data("rnaMetrics")

dataFrame <- stability(data=rnaMetrics, k=2, bs=100, getImages = FALSE)
dataFrame <- stabilityRange(data=rnaMetrics, k.range=c(2,4), bs=20, getImages = FALSE)
dataFrame <- stabilitySet(data=rnaMetrics, k.set=c(2,3,4), bs=20, getImages = FALSE)

dataFrame <- quality(data=rnaMetrics, k=3, getImages = FALSE)
dataFrame <- qualityRange(data=rnaMetrics, k.range=c(2,4), getImages = FALSE)
assay(getDataQualityRange(dataFrame, 2), 1)
dataFrame1 <- qualitySet(data=rnaMetrics, k.set=c(2,3,4), getImages = FALSE)


dataFrame <- metricsCorrelations(data=rnaMetrics, getImages = FALSE, margins = c(4,4,11,10))
assay(dataFrame, 1)

