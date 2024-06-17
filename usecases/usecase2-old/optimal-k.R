library(evaluomeR)

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
inputAgroPath=paste0(wd,"data/agro.csv")
inputOboPath=paste0(wd,"data/obo-119.csv")

agroData = read.csv(inputAgroPath, header = TRUE)
oboData = read.csv(inputOboPath, header = TRUE)
bothData = rbind(agroData, oboData)

agroStab = stabilityRange(data=agroData, k.range=c(2,6), bs=100, getImages = FALSE, seed=13606)
oboStab = stabilityRange(data=oboData, k.range=c(2,6), bs=100, getImages = FALSE, seed=13606)
bothStab = stabilityRange(data=bothData, k.range=c(2,6), bs=100, getImages = FALSE, seed=13606)

agroSil <- qualityRange(data=agroData, k.range=c(2,6), getImages = FALSE, seed=13606)
oboSil <- qualityRange(data=oboData, k.range=c(2,6), getImages = FALSE, seed=13606)
bothSil <- qualityRange(data=bothData, k.range=c(2,6), getImages = FALSE, seed=13606)

agroOpt <- getOptimalKValue(agroStab, agroSil, k.range=c(3,6))
oboOpt <- getOptimalKValue(oboStab, oboSil, k.range=c(3,6))
bothOpt <- getOptimalKValue(bothStab, bothSil, k.range=c(3,6))

agroOpt <- agroOpt[c("Global_optimal_k")]
oboOpt <- oboOpt[c("Global_optimal_k")]
bothOpt <- bothOpt[c("Global_optimal_k")]

kOptTable <- NULL
kOptTable$agro <- agroOpt$Global_optimal_k
kOptTable$obo <- oboOpt$Global_optimal_k
kOptTable$both <- bothOpt$Global_optimal_k
kOptTable <- as.data.frame(kOptTable)
rownames(kOptTable) <- rownames(agroOpt)

View(kOptTable)
