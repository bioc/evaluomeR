library(evaluomeR)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(stringr)
library(reshape2)
library(ggthemes)

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
seed = 13606
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
inputBangertPath=paste0(wd,"data/bangertdrowns2004.csv")
inputLiPath=paste0(wd,"data/li2007.csv")
inputMolloyPath=paste0(wd,"data/molloy2014.csv")

bangertData = read.csv(inputBangertPath, header = TRUE)
liData = read.csv(inputLiPath, header = TRUE)
molloyData = read.csv(inputMolloyPath, header = TRUE)

#### Stability [2,6] ----
k.range=c(2,6)
stabBangert <- stabilityRange(data=bangertData, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabLi <- stabilityRange(data=liData, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabMolloy <- stabilityRange(data=molloyData, k.range=k.range, bs=100, getImages = FALSE, seed=seed)

#### Stability SE to DF ----
meanStabBangert = standardizeStabilityData(stabBangert, k.range)
meanStabBangert$Metric = rownames(meanStabBangert)
meanStabBangert$Repository = "Bangert"
colnames(meanStabBangert) = str_remove_all(colnames(meanStabBangert), "k_")

meanStabLi = standardizeStabilityData(stabLi, k.range)
meanStabLi$Metric = rownames(meanStabLi)
meanStabLi$Repository = "Li"
colnames(meanStabLi) = str_remove_all(colnames(meanStabLi), "k_")

meanStabMolloy = standardizeStabilityData(stabMolloy, k.range)
meanStabMolloy$Metric = rownames(meanStabMolloy)
meanStabMolloy$Repository = "Molloy"
colnames(meanStabMolloy) = str_remove_all(colnames(meanStabMolloy), "k_")

#### Stability plotting ----

stabDfList = list(meanStabBangert, meanStabLi, meanStabMolloy)

stabPlot = ggplot()
min = Inf
max = -Inf
for (stabDf in stabDfList) {
  stabDf_melt = melt(stabDf, id.vars = c("Metric", "Repository"))
  stabDf_melt$variable = as.integer(stabDf_melt$variable)
  stabPlot = stabPlot +
    geom_line(stabDf_melt,
              mapping = aes(x=variable, y=value, group = 1, linetype=Metric, colour=Metric)) +
    geom_point(stabDf_melt,
               mapping = aes(x=variable, y=value, group = 1, shape=Repository))
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

ggsave(plot = stabPlot, filename=paste0(plotDir, "/stability_metafor.pdf"),
       device="pdf", units="cm", width = 20, height = 10, dpi="print")

#### Quality [2,6] ----
qualBangert <- qualityRange(data=bangertData, k.range=k.range, getImages = FALSE, seed=seed)
qualLi <- qualityRange(data=liData, k.range=k.range, getImages = FALSE, seed=seed)
qualMolloy <- qualityRange(data=molloyData, k.range=k.range, getImages = FALSE, seed=seed)

#### Quality SE to DF ----
silBangert = standardizeQualityData(qualBangert, k.range)
silBangert$Metric = rownames(silBangert)
silBangert$Repository = "Bangert"
colnames(silBangert) = str_remove_all(colnames(silBangert), "k_")

silLi = standardizeQualityData(qualLi, k.range)
silLi$Metric = rownames(silLi)
silLi$Repository = "Li"
colnames(silLi) = str_remove_all(colnames(silLi), "k_")

silMolloy = standardizeQualityData(qualMolloy, k.range)
silMolloy$Metric = rownames(silMolloy)
silMolloy$Repository = "Molloy"
colnames(silMolloy) = str_remove_all(colnames(silMolloy), "k_")

#### Quality plotting ----

silDfList = list(silBangert, silLi, silMolloy)

silPlot = ggplot()
min = Inf
max = -Inf
for (silDf in silDfList) {
  silDf_melt = melt(silDf, id.vars = c("Metric", "Repository"))
  silDf_melt$variable = as.integer(silDf_melt$variable)
  silPlot = silPlot +
    geom_line(silDf_melt,
              mapping = aes(x=variable, y=value, group = 1, linetype=Metric, colour=Metric)) +
    geom_point(silDf_melt,
               mapping = aes(x=variable, y=value, group = 1, shape=Repository))
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

ggsave(plot = silPlot, filename=paste0(plotDir, "/silhouette_metafor.pdf"),
       device="pdf", units="cm", width = 20, height = 10, dpi="print")

pg = plot_grid(stabPlot, silPlot, align = "v", nrow = 2, rel_heights = c(1/2, 1/2), labels=c("A", "B"))
save_plot(paste0(plotDir, "/stability_silhouette_metafor.pdf"), pg, nrow=2, dpi="print")
save_plot(paste0(plotDir, "/stability_silhouette_metafor.tiff"), pg, nrow=2, dpi="print", device="tiff")

#### CSV generation ----

# Stab
write.csv(standardizeStabilityData(stabBangert, k.range), paste0(dataDir, "/stabilityBangert.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabLi, k.range), paste0(dataDir, "/stabilityLi.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabMolloy, k.range), paste0(dataDir, "/stabilityMolloy.csv"), row.names = TRUE)

# Qual
write.csv(standardizeQualityData(qualBangert, k.range), paste0(dataDir, "/qualityBangert.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualLi, k.range), paste0(dataDir, "/qualityLi.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualMolloy, k.range), paste0(dataDir, "/qualityMolloy.csv"), row.names = TRUE)

#### Optimal k values ----
k.range=c(3,6)

# Bangert
kOptBangert = getOptimalKValue(stabBangert, qualBangert, k.range=k.range)$Global_optimal_k
# Li
kOptLi = getOptimalKValue(stabLi, qualLi, k.range=k.range)$Global_optimal_k
# Molloy
kOptMolloy = getOptimalKValue(stabMolloy, qualMolloy, k.range=k.range)$Global_optimal_k

tableKOpt = setNames(as.data.frame(matrix(nrow=1,ncol=3)), c("Bangert","Li","Molloy"))
rownames(tableKOpt) <- "Metric"
tableKOpt["Bangert"] <- kOptBangert
tableKOpt["Li"] <- kOptLi
tableKOpt["Molloy"] <- kOptMolloy
tableKOpt

