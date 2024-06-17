library(evaluomeR)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(stringr)
library(reshape2)
library(ggthemes)

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
seed = 13000
plotDir=paste0(wd,"plots")
dataDir=paste0(wd,"results-csv")

dir.create(file.path(plotDir))
dir.create(file.path(dataDir))

getFormattedK <- function(k) {
  return(gsub("^.*_","", k))
}

standardizeStabilityData <- function(stabData, k.range=NULL) {
  stabDf = as.data.frame(assay(stabData)) # Getting first assay, which is 'stabData$stability_mean'
  lengthColnames = length(colnames(stabDf))
  toRemove = list()
  for (i in seq(1, lengthColnames, 1)) {
    colname = colnames(stabDf)[i]
    newColname = gsub("^.*_.*_.*_","k_", colname)
    colnames(stabDf)[i] = newColname
    if (i != 1) { # Skip Metric column
      k = as.numeric(getFormattedK(newColname))
      if (!is.null(k.range) && (k < k.range[1] || k > k.range[2])) {
        toRemove = append(toRemove, newColname)
        next
      }
      stabDf[newColname] = as.numeric(as.character(stabDf[[newColname]]))
    }
  }

  for (columnName in toRemove) {
    stabDf[, columnName] = list(NULL)
    lengthColnames = lengthColnames-1
  }

  inputStartRange = as.numeric(getFormattedK(colnames(stabDf)[2]))
  inputEndRange = as.numeric(getFormattedK(colnames(stabDf)[lengthColnames]))
  if (!is.null(k.range) && (k.range[1] < inputStartRange || k.range[2] > inputEndRange)) {
    # Input k.range is not a subset of the stabData k ranges
    stop("Input k.range [", k.range[1], ", ", k.range[2], "] is not a subset of data range [",
         inputStartRange, ", ", inputEndRange, "]")
  }

  rownames(stabDf) = stabDf$Metric
  stabDf = stabDf[, -1] # Remove "Metric" column, metrics are rownames now
  stabDf <- stabDf[ order(row.names(stabDf)), ]
  return(stabDf)
}

standardizeQualityData <- function(qualData, k.range=NULL) {
  lengthQuality = length(qualData)
  qualRangeStart = getFormattedK(names(qualData)[1])
  qualRangeEnd = getFormattedK(names(qualData)[lengthQuality])
  Metric = NULL
  kValues = list()
  for (i in seq(qualRangeStart, qualRangeEnd, 1)) {
    curQual = as.data.frame(assay(getDataQualityRange(qualData, i)))
    if (i == qualRangeStart) {
      Metric = as.character(curQual$Metric)
    }
    kValues[[i]] = as.numeric(as.character(curQual$Avg_Silhouette_Width))
  }
  qualDf = as.data.frame(Metric)
  for (i in seq(qualRangeStart, qualRangeEnd, 1)) {
    values = kValues[[i]]
    newColname = paste0("k_", i)
    k = as.numeric(getFormattedK(newColname))
    if (!is.null(k.range) && (k < k.range[1] || k > k.range[2])) {
      next
    }
    qualDf[[newColname]] = values
  }

  if (!is.null(k.range) && (k.range[1] < qualRangeStart || k.range[2] > qualRangeEnd)) {
    # Input k.range is not a subset of the stabData k ranges
    stop("Input k.range [", k.range[1], ", ", k.range[2], "] is not a subset of range [",
         qualRangeStart, ", ", qualRangeEnd, "]")
  }

  rownames(qualDf) = qualDf$Metric
  qualDf = qualDf[, -1] # Remove "Metric" column, metrics are rownames now
  qualDf <- qualDf[ order(row.names(qualDf)), ]
  return(qualDf)
}

#### Data loading ----
matcomp15Path=paste0(wd,"data/data-matcomp15.csv")
matcomp16Path=paste0(wd,"data/data-matcomp16.csv")
matcomp17Path=paste0(wd,"data/data-matcomp17.csv")

multi15Path=paste0(wd,"data/data-multi15.csv")
multi16Path=paste0(wd,"data/data-multi16.csv")
multi17Path=paste0(wd,"data/data-multi17.csv")

multimc15Path=paste0(wd,"data/data-multimc15.csv")
multimc16Path=paste0(wd,"data/data-multimc16.csv")
multimc17Path=paste0(wd,"data/data-multimc17.csv")

inputMatcomp15 = read.csv(matcomp15Path, header = TRUE)
inputMatcomp16 = read.csv(matcomp16Path, header = TRUE)
inputMatcomp17 = read.csv(matcomp17Path, header = TRUE)

inputMulti15 = read.csv(multi15Path, header = TRUE)
inputMulti16 = read.csv(multi16Path, header = TRUE)
inputMulti17 = read.csv(multi17Path, header = TRUE)

inputMultimc15 = read.csv(multimc15Path, header = TRUE)
inputMultimc16 = read.csv(multimc16Path, header = TRUE)
inputMultimc17 = read.csv(multimc17Path, header = TRUE)

#### Stability [2,15] ----
k.range=c(2,15)
stabMatcomp15 <- stabilityRange(data=inputMatcomp15,
                                k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabMatcomp16 <- stabilityRange(data=inputMatcomp16, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabMatcomp17 <- stabilityRange(data=inputMatcomp17, k.range=k.range, bs=100, getImages = FALSE, seed=seed)

stabMulti15 <- stabilityRange(data=inputMulti15, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabMulti16 <- stabilityRange(data=inputMulti16, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabMulti17 <- stabilityRange(data=inputMultimc17, k.range=k.range, bs=100, getImages = FALSE, seed=seed)

stabMultimc15 <- stabilityRange(data=inputMultimc15, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabMultimc16 <- stabilityRange(data=inputMultimc16, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabMultimc17 <- stabilityRange(data=inputMultimc17, k.range=k.range, bs=100, getImages = FALSE, seed=seed)

#### Stability SE to DF ----
meanStabMatcomp15 = standardizeStabilityData(stabMatcomp15, k.range)
meanStabMatcomp15$Metric = "Matcomp15"
meanStabMatcomp15$category = "MCB"
meanStabMatcomp15$year = "2015"
colnames(meanStabMatcomp15) = str_remove_all(colnames(meanStabMatcomp15), "k_")

meanStabMatcomp16 = standardizeStabilityData(stabMatcomp16, k.range)
meanStabMatcomp16$Metric = "Matcomp16"
meanStabMatcomp16$category = "MCB"
meanStabMatcomp16$year = "2016"
colnames(meanStabMatcomp16) = str_remove_all(colnames(meanStabMatcomp16), "k_")

meanStabMatcomp17 = standardizeStabilityData(stabMatcomp17, k.range)
meanStabMatcomp17$Metric = "Matcomp17"
meanStabMatcomp17$category = "MCB"
meanStabMatcomp17$year = "2017"
colnames(meanStabMatcomp17) = str_remove_all(colnames(meanStabMatcomp17), "k_")

meanStabMulti15 = standardizeStabilityData(stabMulti15, k.range)
meanStabMulti15$Metric = "Multi15"
meanStabMulti15$category = "MS"
meanStabMulti15$year = "2015"
colnames(meanStabMulti15) = str_remove_all(colnames(meanStabMulti15), "k_")

meanStabMulti16 = standardizeStabilityData(stabMulti16, k.range)
meanStabMulti16$Metric = "Multi16"
meanStabMulti16$category = "MS"
meanStabMulti16$year = "2016"
colnames(meanStabMulti16) = str_remove_all(colnames(meanStabMulti16), "k_")

meanStabMulti17 = standardizeStabilityData(stabMulti17, k.range)
meanStabMulti17$Metric = "Multi17"
meanStabMulti17$category = "MS"
meanStabMulti17$year = "2017"
colnames(meanStabMulti17) = str_remove_all(colnames(meanStabMulti17), "k_")

meanStabMultimc15 = standardizeStabilityData(stabMultimc15, k.range)
meanStabMultimc15$Metric = "Multimc15"
meanStabMultimc15$category = "MCB+MS"
meanStabMultimc15$year = "2015"
colnames(meanStabMultimc15) = str_remove_all(colnames(meanStabMultimc15), "k_")

meanStabMultimc16 = standardizeStabilityData(stabMultimc16, k.range)
meanStabMultimc16$Metric = "Multimc16"
meanStabMultimc16$category = "MCB+MS"
meanStabMultimc16$year = "2016"
colnames(meanStabMultimc16) = str_remove_all(colnames(meanStabMultimc16), "k_")

meanStabMultimc17 = standardizeStabilityData(stabMultimc17, k.range)
meanStabMultimc17$Metric = "Multimc17"
meanStabMultimc17$category = "MCB+MS"
meanStabMultimc17$year = "2017"
colnames(meanStabMultimc17) = str_remove_all(colnames(meanStabMultimc17), "k_")

#### Stability plotting ----

stabDfList = list(meanStabMatcomp15, meanStabMatcomp16, meanStabMatcomp17,
                  meanStabMulti15, meanStabMulti16, meanStabMulti17,
                  meanStabMultimc15, meanStabMultimc16, meanStabMultimc17)

stabPlot = ggplot()
min = Inf
max = -Inf
for (stabDf in stabDfList) {
  stabDf_melt = melt(stabDf, id.vars = c("Metric", "category", "year"))
  colnames(stabDf_melt)[1] = "Journal"
  stabDf_melt$variable = as.integer(stabDf_melt$variable)
  stabPlot = stabPlot +
    geom_line(stabDf_melt,
              mapping = aes(x=variable, y=value, group = 1, linetype=category, colour=category)) +
    geom_point(stabDf_melt,
               mapping = aes(x=variable, y=value, group = 1, shape=year))
  if (min > min(stabDf_melt$value)) {
    min = min(stabDf_melt$value)
  }
  if (max < max(stabDf_melt$value)) {
    max = max(stabDf_melt$value)
  }
}

stabPlot = stabPlot +
  scale_x_continuous(name="k", breaks=1:14, labels=2:15) +
  scale_y_continuous(name="Stability", limits = c(min,max), breaks = c(0.6, 0.75, 0.85), labels=c(0.6, 0.75, 0.85)) +
  scale_colour_grey(start = 0.7, end = 0) +
  theme_calc(base_family = "sans")

ggsave(plot = stabPlot, filename=paste0(plotDir, "/stability_impact_factor.pdf"),
       device="pdf", units="cm", width = 20, height = 10, dpi="print")

#### Quality [2,15] ----
qualMatcomp15 <- qualityRange(data=inputMatcomp15, k.range=k.range, getImages = FALSE, seed=seed)
qualMatcomp16 <- qualityRange(data=inputMatcomp16, k.range=k.range, getImages = FALSE, seed=seed)
qualMatcomp17 <- qualityRange(data=inputMatcomp17, k.range=k.range, getImages = FALSE, seed=seed)

qualMulti15 <- qualityRange(data=inputMulti15, k.range=k.range, getImages = FALSE, seed=seed)
qualMulti16 <- qualityRange(data=inputMulti16, k.range=k.range, getImages = FALSE, seed=seed)
qualMulti17 <- qualityRange(data=inputMultimc17, k.range=k.range, getImages = FALSE, seed=seed)

qualMultimc15 <- qualityRange(data=inputMultimc15, k.range=k.range, getImages = FALSE, seed=seed)
qualMultimc16 <- qualityRange(data=inputMultimc16, k.range=k.range, getImages = FALSE, seed=seed)
qualMultimc17 <- qualityRange(data=inputMultimc17, k.range=k.range, getImages = FALSE, seed=seed)

#### Quality SE to DF ----
silMatcomp15 = standardizeQualityData(qualMatcomp15, k.range)
silMatcomp15$Metric = "Matcomp15"
silMatcomp15$category = "MCB"
silMatcomp15$year = "2015"
colnames(silMatcomp15) = str_remove_all(colnames(silMatcomp15), "k_")

silMatcomp16 = standardizeQualityData(qualMatcomp16, k.range)
silMatcomp16$Metric = "Matcomp16"
silMatcomp16$category = "MCB"
silMatcomp16$year = "2016"
colnames(silMatcomp16) = str_remove_all(colnames(silMatcomp16), "k_")

silMatcomp17 = standardizeQualityData(qualMatcomp17, k.range)
silMatcomp17$Metric = "Matcomp17"
silMatcomp17$category = "MCB"
silMatcomp17$year = "2017"
colnames(silMatcomp17) = str_remove_all(colnames(silMatcomp17), "k_")

silMulti15 = standardizeQualityData(qualMulti15, k.range)
silMulti15$Metric = "Multi15"
silMulti15$category = "MS"
silMulti15$year = "2015"
colnames(silMulti15) = str_remove_all(colnames(silMulti15), "k_")

silMulti16 = standardizeQualityData(qualMulti16, k.range)
silMulti16$Metric = "Multi16"
silMulti16$category = "MS"
silMulti16$year = "2016"
colnames(silMulti16) = str_remove_all(colnames(silMulti16), "k_")

silMulti17 = standardizeQualityData(qualMulti17, k.range)
silMulti17$Metric = "Multi17"
silMulti17$category = "MS"
silMulti17$year = "2017"
colnames(silMulti17) = str_remove_all(colnames(silMulti17), "k_")

silMultimc15 = standardizeQualityData(qualMultimc15, k.range)
silMultimc15$Metric = "Multimc15"
silMultimc15$category = "MCB+MS"
silMultimc15$year = "2015"
colnames(silMultimc15) = str_remove_all(colnames(silMultimc15), "k_")

silMultimc16 = standardizeQualityData(qualMultimc16, k.range)
silMultimc16$Metric = "Multimc16"
silMultimc16$category = "MCB+MS"
silMultimc16$year = "2016"
colnames(silMultimc16) = str_remove_all(colnames(silMultimc16), "k_")

silMultimc17 = standardizeQualityData(qualMultimc16, k.range)
silMultimc17$Metric = "Multimc17"
silMultimc17$category = "MCB+MS"
silMultimc17$year = "2017"
colnames(silMultimc17) = str_remove_all(colnames(silMultimc17), "k_")

#### Quality plotting ----

silDfList = list(silMatcomp15, silMatcomp16, silMatcomp17,
                  silMulti15, silMulti16, silMulti17,
                  silMultimc15, silMultimc16, silMultimc17)

silPlot = ggplot()
min = Inf
max = -Inf
for (silDf in silDfList) {
  silDf_melt = melt(silDf, id.vars = c("Metric", "category", "year"))
  colnames(silDf_melt)[1] = "Journal"
  silDf_melt$variable = as.integer(silDf_melt$variable)
  silPlot = silPlot +
    geom_line(silDf_melt,
              mapping = aes(x=variable, y=value, group = 1, linetype=category, colour=category)) +
    geom_point(silDf_melt,
               mapping = aes(x=variable, y=value, group = 1, shape=year))
  if (min > min(silDf_melt$value)) {
    min = min(silDf_melt$value)
  }
  if (max < max(silDf_melt$value)) {
    max = max(silDf_melt$value)
  }
}

silPlot = silPlot +
  scale_x_continuous(name="k", breaks=1:14, labels=2:15) +
  scale_y_continuous(name="Quality", limits = c(min,max), breaks = c(0.25, 0.5, 0.7), labels=c(0.25, 0.5, 0.7)) +
  scale_colour_grey(start = 0.7, end = 0) +
  theme_calc()

ggsave(plot = silPlot, filename=paste0(plotDir, "/silhouette_impact_factor.pdf"),
       device="pdf", units="cm", width = 20, height = 10, dpi="print")

pg = plot_grid(stabPlot, silPlot, align = "v", nrow = 2, rel_heights = c(1/2, 1/2), labels=c("A", "B"))
save_plot(paste0(plotDir, "/stability_silhouette_impact_factor.pdf"), pg, nrow=2, dpi="print")
save_plot(paste0(plotDir, "/stability_silhouette_impact_factor.tiff"), pg, nrow=2, dpi="print", device="tiff")

#### CSV generation ----

# MC Stab
write.csv(standardizeStabilityData(stabMatcomp15, k.range), paste0(dataDir, "/stabilityMatcomp15.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabMatcomp16, k.range), paste0(dataDir, "/stabilityMatcomp16.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabMatcomp17, k.range), paste0(dataDir, "/stabilityMatcomp17.csv"), row.names = TRUE)

# MC Qual
write.csv(standardizeQualityData(qualMatcomp15, k.range), paste0(dataDir, "/qualityMatcomp15.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualMatcomp16, k.range), paste0(dataDir, "/qualityMatcomp16.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualMatcomp17, k.range), paste0(dataDir, "/qualityMatcomp17.csv"), row.names = TRUE)

# Multi Stab
write.csv(standardizeStabilityData(stabMulti15, k.range), paste0(dataDir, "/stabilityMulti15.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabMulti16, k.range), paste0(dataDir, "/stabilityMulti16.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabMulti17, k.range), paste0(dataDir, "/stabilityMulti17.csv"), row.names = TRUE)

# Multi Qual
write.csv(standardizeQualityData(qualMulti15, k.range), paste0(dataDir, "/qualityMulti15.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualMulti16, k.range), paste0(dataDir, "/qualityMulti16.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualMulti17, k.range), paste0(dataDir, "/qualityMulti17.csv"), row.names = TRUE)

# Multi + MC Stab
write.csv(standardizeStabilityData(stabMultimc15, k.range), paste0(dataDir, "/stabilityMultimc15.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabMultimc16, k.range), paste0(dataDir, "/stabilityMultimc16.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabMultimc17, k.range), paste0(dataDir, "/stabilityMultimc17.csv"), row.names = TRUE)

# Multi + MC Qual
write.csv(standardizeQualityData(qualMultimc15, k.range), paste0(dataDir, "/qualityMultimc15.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualMultimc16, k.range), paste0(dataDir, "/qualityMultimc16.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualMultimc17, k.range), paste0(dataDir, "/qualityMultimc17.csv"), row.names = TRUE)

#### Optimal k values ----
k.range=c(3,15)

# MC
kOptMC15 = getOptimalKValue(stabMatcomp15, qualMatcomp15, k.range=k.range)$Global_optimal_k
kOptMC16 = getOptimalKValue(stabMatcomp16, qualMatcomp16, k.range=k.range)$Global_optimal_k
kOptMC17 = getOptimalKValue(stabMatcomp17, qualMatcomp17, k.range=k.range)$Global_optimal_k

# Multi
kOptMulti15 = getOptimalKValue(stabMulti15, qualMulti15, k.range=k.range)$Global_optimal_k
kOptMulti16 = getOptimalKValue(stabMulti16, qualMulti16, k.range=k.range)$Global_optimal_k
kOptMulti17 = getOptimalKValue(stabMulti17, qualMulti17, k.range=k.range)$Global_optimal_k

# Multi + MC
kOptMultimc15 = getOptimalKValue(stabMultimc15, qualMultimc15, k.range=k.range)$Global_optimal_k
kOptMultimc16 = getOptimalKValue(stabMultimc16, qualMultimc16, k.range=k.range)$Global_optimal_k
kOptMultimc17 = getOptimalKValue(stabMultimc17, qualMultimc17, k.range=k.range)$Global_optimal_k

tableKOpt = setNames(as.data.frame(matrix(nrow=3,ncol=3)), c("2015","2016","2017"))
rownames(tableKOpt) <- c("MC", "Multi", "Multimc")
tableKOpt["2015"] <- c(kOptMC15, kOptMulti15, kOptMultimc15)
tableKOpt["2016"] <- c(kOptMC16, kOptMulti16, kOptMultimc16)
tableKOpt["2017"] <- c(kOptMC17, kOptMulti17, kOptMultimc17)
tableKOpt
