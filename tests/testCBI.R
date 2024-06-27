library(evaluomeR)


evaluomeRSupportedCBI()


dataFrame <- stability(data=ontMetrics, cbi="kmeans", k=3, all_metrics=FALSE, bs=100)
assay(dataFrame)

#dataFrame <- stabilityRange(data=ontMetrics, cbi="rskc", k.range=c(3,4), all_metrics=TRUE, bs=100, L1=2)
#assay(dataFrame)

#dataFrame <- stabilitySet(data=ontMetrics, k.set=c(3,4), bs=100, cbi="rskc", all_metrics=TRUE, L1=2)
#assay(dataFrame)

#dataFrame <- quality(data=ontMetrics, cbi="rskc", k=3, all_metrics=TRUE, L1=2)
#assay(dataFrame)

#dataFrame <- qualityRange(data=ontMetrics, cbi="rskc", k.range=c(3,4), all_metrics=TRUE, L1=2)
#assay(dataFrame$k_3)

#dataFrame <- qualitySet(data=ontMetrics, cbi="rskc", k.set=c(3,5), all_metrics=TRUE, L1=2)
#assay(dataFrame$k_3)

