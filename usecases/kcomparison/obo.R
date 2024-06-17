library(evaluomeR)
source(paste0(wd,"functions.R"))

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
inputPath=paste0(wd,"data/obo-119.csv")
outputDir=paste0(wd,"results-obo")

inputDataObo = read.csv(inputPath, header = TRUE)

correlationsOBO = metricsCorrelations(inputDataObo)
correlationsOBO = assay(correlationsOBO)

plotMetricsMinMax(inputDataObo) # minmax
plotMetricsBoxplot(inputDataObo) # boxplot
plotMetricsCluster(inputDataObo) # cluster
plotMetricsViolin(inputDataObo) # violin

if (!exists("GLOBAL_CBI")) {
  GLOBAL_CBI = "kmeans"
}

bs=500

cat("Using CBI: ", GLOBAL_CBI, "\n")

stabilityDataObo <- stabilityRange(data=inputDataObo,
                                   k.range=c(3,15), cbi=GLOBAL_CBI, bs=bs,
                                   getImages = FALSE, seed=13606)
qualityDataObo <- qualityRange(data=inputDataObo,
                               k.range=c(3,15), cbi=GLOBAL_CBI, getImages = FALSE,
                               seed=13606)

kOptTableOboFull <- getOptimalKValue(stabilityDataObo, qualityDataObo, k.range=c(3,15))
kOptTableObo <- kOptTableOboFull[, c("Stability_max_k", "Quality_max_k", "Global_optimal_k")]
rownames(kOptTableObo) <- kOptTableOboFull$Metric
kOptTableOboValues <- kOptTableOboFull[, c("Stability_max_k_stab", "Quality_max_k_stab")]
csvPath = paste0(outputDir, "/clusters/obo.optimal.k_",GLOBAL_CBI,".csv")
#write.csv(kOptTableObo, csvPath, row.names = TRUE)
#> write.csv(kOptTableOboFull, csvPath, row.names = TRUE)
#View(kOptTableObo)

stabilityTableObo <- getDfStabilityData(stabilityDataObo)
silhouetteTableObo <- getDfQualityData(qualityDataObo)

csvPath = paste0(outputDir, "/clusters/obo.stability_",GLOBAL_CBI,".csv")
#write.csv(stabilityTableObo, csvPath, row.names = TRUE)
csvPath = paste0(outputDir, "/clusters/obo.goodness_",GLOBAL_CBI,".csv")
#write.csv(silhouetteTableObo, csvPath, row.names = TRUE)

# Conteo de scores estables:

length(which(stabilityTableObo["k_3"]>0.75)) # At least stable, (0.75,0.85]
length(which(silhouetteTableObo["k_3"]>0.5)) # At least reasonable structure, (0.50, 0.70]
