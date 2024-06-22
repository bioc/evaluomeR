library(evaluomeR)

data("rnaMetrics")

dataFrame <- stability(data=rnaMetrics, k=4, bs=100, all_metrics = FALSE, getImages = FALSE)
dataFrame <- stabilityRange(data=rnaMetrics, k.range=c(2,4), bs=20, all_metrics = FALSE, getImages = FALSE)
assay(dataFrame)
# Metric    Mean_stability_k_2  Mean_stability_k_3  Mean_stability_k_4
# [1,] "RIN"     "0.825833333333333" "0.778412698412698" "0.69625"
# [2,] "DegFact" "0.955595238095238" "0.977777777777778" "0.820833333333333"
dataFrame <- stabilitySet(data=rnaMetrics, k.set=c(2,3,4), bs=20, all_metrics = FALSE, getImages = FALSE)

dataFrame <- quality(data=rnaMetrics, cbi="kmeans", k=3, all_metrics = FALSE, getImages = FALSE)
assay(dataFrame)
# Metric    Cluster_1_SilScore  Cluster_2_SilScore  Cluster_3_SilScore
# [1,] "RIN"     "0.420502645502646" "0.724044583696066" "0.68338517747747"
# [2,] "DegFact" "0.876516605981734" "0.643613928123002" "0.521618857725795"
# Avg_Silhouette_Width Cluster_1_Size Cluster_2_Size Cluster_3_Size
# [1,] "0.627829396038413"  "4"            "4"            "8"
# [2,] "0.737191191352892"  "8"            "5"            "3"
dataFrame <- qualityRange(data=rnaMetrics, k.range=c(2,4), seed = 20, all_metrics = FALSE, getImages = FALSE)
assay(getDataQualityRange(dataFrame, 2))
# Metric    Cluster_1_SilScore  Cluster_2_SilScore  Avg_Silhouette_Width Cluster_1_Size
# 1 "RIN"     "0.583166775069983" "0.619872562681118" "0.608402004052639"  "5"
# 2 "DegFact" "0.664573423022171" "0.675315791048653" "0.666587617027136"  "13"
# Cluster_2_Size
# 1 "11"
# 2 "3"
assay(getDataQualityRange(dataFrame, 4))
# Metric    Cluster_1_SilScore  Cluster_2_SilScore  Cluster_3_SilScore
# 1 "RIN"     "0.420502645502646" "0.674226581940152" "0.433333333333333"
# 2 "DegFact" "0.759196481622952" "0.59496499852177"  "0.600198799385732"
# Cluster_4_SilScore  Avg_Silhouette_Width Cluster_1_Size Cluster_2_Size Cluster_3_Size
# 1 "0.348714574898785" "0.463905611516569"  "4"            "4"            "3"
# 2 "0.521618857725795" "0.634170498361632"  "5"            "3"            "5"
# Cluster_4_Size
# 1 "5"
# 2 "3"
dataFrame1 <- qualitySet(data=rnaMetrics, k.set=c(2,3,4), all_metrics = FALSE, getImages = FALSE)


dataFrame <- metricsCorrelations(data=rnaMetrics, getImages = FALSE, margins = c(4,4,11,10))
assay(dataFrame, 1)


dataFrame <- stability(data=rnaMetrics, cbi="kmeans", k=2, bs=100, all_metrics = FALSE, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="clara", k=2, bs=100, all_metrics = FALSE, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="clara_pam", k=2, bs=100, all_metrics = FALSE, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="hclust", k=2, bs=100, all_metrics = FALSE, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="pamk", k=2, bs=100, all_metrics = FALSE, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="pamk_pam", k=2, bs=100, all_metrics = FALSE, getImages = FALSE)
#dataFrame <- stability(data=rnaMetrics, cbi="rskc", k=2, bs=100, all_metrics = TRUE, L1 = 2, alpha=0, getImages = FALSE)

# Supported CBIs:
evaluomeRSupportedCBI()

dataFrame <- qualityRange(data=rnaMetrics, k.range=c(2,10), all_metrics = FALSE, getImages = FALSE)
dataFrame

#dataFrame <- stabilityRange(data=rnaMetrics, k.range=c(2,8), bs=20, getImages = FALSE)
#assay(dataFrame)

