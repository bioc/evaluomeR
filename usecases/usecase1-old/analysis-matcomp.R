library(evaluomeR)

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
input15Path=paste0(wd,"data/data-matcomp15.csv")
input16Path=paste0(wd,"data/data-matcomp16.csv")
input17Path=paste0(wd,"data/data-matcomp17.csv")
outputDir=paste0(wd,"results-matcomp")

inputData15 = read.csv(input15Path, header = TRUE)
inputData16 = read.csv(input16Path, header = TRUE)
inputData17 = read.csv(input17Path, header = TRUE)

# Stability [2,15]
dataFrame15 <- stabilityRange(data=inputData15, k.range=c(2,15), bs=100, getImages = TRUE)
dataFrame16 <- stabilityRange(data=inputData16, k.range=c(2,15), bs=100, getImages = TRUE)
dataFrame17 <- stabilityRange(data=inputData17, k.range=c(2,15), bs=100, getImages = TRUE)
dataStability15 = assay(dataFrame15, 1)
dataStability16 = assay(dataFrame16, 1)
dataStability17 = assay(dataFrame17, 1)
csvPath = paste0(outputDir, "/stabilitymc15.csv")
write.csv(dataStability15, csvPath, row.names = FALSE)
csvPath = paste0(outputDir, "/stabilitymc16.csv")
write.csv(dataStability16, csvPath, row.names = FALSE)
csvPath = paste0(outputDir, "/stabilitymc17.csv")
write.csv(dataStability17, csvPath, row.names = FALSE)

# Quality [2,15]
dataFrame15 <- qualityRange(data=inputData15, k.range=c(2,15), getImages = TRUE)
dataFrame16 <- qualityRange(data=inputData16, k.range=c(2,15), getImages = TRUE)
dataFrame17 <- qualityRange(data=inputData17, k.range=c(2,15), getImages = TRUE)
for (i in 2:15) {
  data15 = assay(getDataQualityRange(dataFrame15, i), 1)
  data16 = assay(getDataQualityRange(dataFrame16, i), 1)
  data17 = assay(getDataQualityRange(dataFrame17, i), 1)
  csvPath = paste0(outputDir, "/qualitymc15_k=", i,".csv")
  write.csv(data15, csvPath, row.names = FALSE)
  csvPath = paste0(outputDir, "/qualitymc16_k=", i,".csv")
  write.csv(data16, csvPath, row.names = FALSE)
  csvPath = paste0(outputDir, "/qualitymc17_k=", i,".csv")
  write.csv(data17, csvPath, row.names = FALSE)
}
