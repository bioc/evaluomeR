library(evaluomeR)

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
source(paste0(wd,"functions.R"))
inputPath=paste0(wd,"data/agro.csv")
outputDir=paste0(wd,"results-agro")

inputDataAgro = read.csv(inputPath, header = TRUE)

correlationsAgro = metricsCorrelations(inputDataAgro)
correlationsAgro = assay(correlationsAgro)

plotMetricsMinMax(inputDataAgro) # minmax
plotMetricsBoxplot(inputDataAgro) # boxplot
plotMetricsCluster(inputDataAgro) # cluster
plotMetricsViolin(inputDataAgro) # violin

#stab = stabilityRange(data=rnaMetrics, k.range=c(3,4), bs=100, getImages=F, seed=13606)
#qual = qualityRange(data=rnaMetrics, k.range=c(3,4), getImages=F, seed=13606)

if (!exists("GLOBAL_CBI")) {
  GLOBAL_CBI = "kmeans"
}

bs=500

cat("Using CBI: ", GLOBAL_CBI, "\n")


# Stability for range [3,8]
stabilityDataAgro <- stabilityRange(data=inputDataAgro,
                                    cbi=GLOBAL_CBI,
                                    k.range=c(3,8), bs=bs,
                                    getImages = FALSE, seed=13606)
# Stability for range [9, 15], since some metrics will be stuck in the
# bootstrap due to a lack of resampling in their values.
newInputDataAgro = inputDataAgro
metricsToRemove = c("AROnto", "DITOnto", "TMOnto2")
newInputDataAgro[, metricsToRemove] <- list(NULL)
newStabilityDataAgro <- stabilityRange(data=newInputDataAgro,
                                    k.range=c(9,15), bs=bs,
                                    getImages = FALSE, seed=13606)

# Stability for DITonto in range [9,14]
ditOntoInputDataAgro = inputDataAgro[, c("Description", "DITOnto")]
ditOntoStabilityDataAgro <- stabilityRange(data=ditOntoInputDataAgro,
                                       k.range=c(9,14), bs=bs,
                                       getImages = FALSE, seed=13606)
# Merge both tables into one
stabilityRange3_15 = getDfStabilityData(stabilityDataAgro, digits = 10)
stabilityRange9_15 = getDfStabilityData(newStabilityDataAgro, digits = 10)
ditOntoStabilityRange9_14 = getDfStabilityData(ditOntoStabilityDataAgro, digits = 10)
missingRowsLength = length(9:15)
for (metric in row.names(stabilityRange3_15)) {
  if (metric == "DITOnto") {
    noNaValues = stabilityRange3_15["DITOnto", ]
    noNaValues = noNaValues[!is.na(noNaValues)]
    stabilityRange3_15[metric, ] = unlist(c(noNaValues, ditOntoStabilityRange9_14["DITOnto", ], NA))
    print(stabilityRange3_15[metric, ])
  } else {
    for (k in 9:15) {
      columnName = paste0("k_",k)
      if (metric %in% metricsToRemove) {
        stabilityRange3_15[metric, ] = unlist(c(stabilityRange3_15[metric, ], rep(NA, missingRowsLength)))
      } else {
        stabilityRange3_15[metric, columnName] = stabilityRange9_15[metric, columnName]
      }
    }
  }
}

stabilityRange3_15NoNa = stabilityRange3_15 # To replace NA's to 0
stabilityRange3_15NoNa[is.na(stabilityRange3_15NoNa)] <- 0
stabilityRange3_15NoNa = as.data.frame(append(stabilityRange3_15NoNa, list(Metric = rownames(stabilityRange3_15)), after = 0))
rownames(stabilityRange3_15NoNa) <- stabilityRange3_15$Metric
colnames(stabilityRange3_15NoNa) <- c("Metric", "Mean_stability_k_3", "Mean_stability_k_4", "Mean_stability_k_5", "Mean_stability_k_6",
                 "Mean_stability_k_7", "Mean_stability_k_8", "Mean_stability_k_9", "Mean_stability_k_10",
                 "Mean_stability_k_11", "Mean_stability_k_12", "Mean_stability_k_13", "Mean_stability_k_14",
                 "Mean_stability_k_15")

qualityDataAgro <- qualityRange(data=inputDataAgro,
                                cbi=GLOBAL_CBI,
                                k.range=c(3,15), getImages = FALSE,
                                seed=13606)


kOptTableAgroFull <- getOptimalKValue(stabilityRange3_15NoNa, qualityDataAgro, k.range=c(3,15))
kOptTableAgro <- kOptTableAgroFull[, c("Stability_max_k", "Quality_max_k", "Global_optimal_k")]
kOptTableAgroValues <- kOptTableAgroFull[, c("Stability_max_k_stab", "Quality_max_k_stab")]
csvPath = paste0(outputDir, "/clusters/agro.optimal.k_",GLOBAL_CBI,".csv")
#write.csv(kOptTableAgro, csvPath, row.names = TRUE)
#> write.csv(kOptTableAgroFull, csvPath, row.names = TRUE, quote = TRUE)
#View(kOptTableAgro)


stabilityTableAgro <- getDfStabilityData(stabilityRange3_15NoNa)
silhouetteTableAgro <- getDfQualityData(qualityDataAgro)

csvPath = paste0(outputDir, "/clusters/agro.stability_",GLOBAL_CBI,".csv")
#write.csv(stabilityTableAgro, csvPath, row.names = TRUE)
csvPath = paste0(outputDir, "/clusters/agro.goodness_",GLOBAL_CBI,".csv")
#write.csv(silhouetteTableAgro, csvPath, row.names = TRUE)

# Conteo de scores estables:

length(which(stabilityTableAgro["k_3"]>0.75)) # At least stable, (0.75,0.85]
length(which(silhouetteTableAgro["k_3"]>0.5)) # At least reasonable structure, (0.50, 0.70]
