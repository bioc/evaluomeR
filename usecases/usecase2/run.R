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
inputAgroPath=paste0(wd,"data/agro.csv")
inputOboPath=paste0(wd,"data/obo-119.csv")

agroData = read.csv(inputAgroPath, header = TRUE)
oboData = read.csv(inputOboPath, header = TRUE)
bothData = rbind(agroData, oboData)


#### Stability [2,6] ----
k.range=c(2,6)
stabAgro <- stabilityRange(data=agroData, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabObo <- stabilityRange(data=oboData, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabBoth <- stabilityRange(data=bothData, k.range=k.range, bs=100, getImages = FALSE, seed=seed)

#### Stability SE to DF ----
meanStabAgro = standardizeStabilityData(stabAgro, k.range)
meanStabAgro$Metric = rownames(meanStabAgro)
meanStabAgro$Repository = "AGRO"
colnames(meanStabAgro) = str_remove_all(colnames(meanStabAgro), "k_")

meanStabObo = standardizeStabilityData(stabObo, k.range)
meanStabObo$Metric = rownames(meanStabObo)
meanStabObo$Repository = "OBO"
colnames(meanStabObo) = str_remove_all(colnames(meanStabObo), "k_")

meanStabBoth = standardizeStabilityData(stabBoth, k.range)
meanStabBoth$Metric = rownames(meanStabBoth)
meanStabBoth$Repository = "AGRO+OBO"
colnames(meanStabBoth) = str_remove_all(colnames(meanStabBoth), "k_")

#### Stability plotting ----
stabDfList = list(meanStabAgro, meanStabObo, meanStabBoth)

stabPlot = ggplot()
min = Inf
max = -Inf
for (stabDf in stabDfList) {
  stabDf = stabDf[which(stabDf$Metric == "ANOnto" |
                          stabDf$Metric == "CBOOnto" |
                          stabDf$Metric == "NOCOnto"
                          ), ]
  stabDf_melt = melt(stabDf, id.vars = c("Metric", "Repository"))
  stabDf_melt$variable = as.integer(stabDf_melt$variable)
  for (metric in c("ANOnto", "CBOOnto", "NOCOnto")) {
    current_melt = stabDf_melt[which(stabDf_melt$Metric == metric), ]
    stabPlot = stabPlot +
      geom_line(current_melt,
                mapping = aes(x=variable, y=value, group = 1, linetype=Metric, colour=Metric)) +
      geom_point(current_melt,
                mapping = aes(x=variable, y=value, group = 1, shape=Repository))
  }

  if (min > min(stabDf_melt$value)) {
    min = min(stabDf_melt$value)
  }
  if (max < max(stabDf_melt$value)) {
    max = max(stabDf_melt$value)
  }
}

stabPlot = stabPlot +
  scale_x_continuous(name="k", breaks=1:5, labels=2:6) +
  scale_y_continuous(name="Stability", limits = c(min,max), breaks = c(0.6, 0.75, 0.85), labels=c(0.6, 0.75, 0.85)) +
  scale_colour_grey(start = 0.7, end = 0) +
  theme_calc()

#ggsave(plot = stabPlot, filename=paste0(plotDir, "/stability_agro_obo_both.pdf"),
#       device="pdf", units="cm", width = 20, height = 10, dpi="print")

#### Quality [2,6] ----
qualAgro <- qualityRange(data=agroData, k.range=k.range, getImages = FALSE, seed=seed)
qualObo <- qualityRange(data=oboData, k.range=k.range, getImages = FALSE, seed=seed)
qualBoth <- qualityRange(data=bothData, k.range=k.range, getImages = FALSE, seed=seed)

#### Quality SE to DF ----
silAgro = standardizeQualityData(qualAgro, k.range)
silAgro$Metric = rownames(silAgro)
silAgro$Repository = "AGRO"
colnames(silAgro) = str_remove_all(colnames(silAgro), "k_")

silObo = standardizeQualityData(qualObo, k.range)
silObo$Metric = rownames(silObo)
silObo$Repository = "OBO"
colnames(silObo) = str_remove_all(colnames(silObo), "k_")

silBoth = standardizeQualityData(qualBoth, k.range)
silBoth$Metric = rownames(silBoth)
silBoth$Repository = "AGRO+OBO"
colnames(silBoth) = str_remove_all(colnames(silBoth), "k_")

#### Quality plotting ----

silDfList = list(silAgro, silObo, silBoth)


silPlot = ggplot()
min = Inf
max = -Inf
for (silDf in silDfList) {
  silDf = silDf[which(silDf$Metric == "ANOnto" |
                        silDf$Metric == "CBOOnto" |
                        silDf$Metric == "NOCOnto"
  ), ]
  silDf_melt = melt(silDf, id.vars = c("Metric", "Repository"))
  silDf_melt$variable = as.integer(silDf_melt$variable)

  for (metric in c("ANOnto", "CBOOnto", "NOCOnto")) {
    current_melt = silDf_melt[which(silDf_melt$Metric == metric), ]
    silPlot = silPlot +
      geom_line(current_melt,
                mapping = aes(x=variable, y=value, group = 1, linetype=Metric, colour=Metric)) +
      geom_point(current_melt,
                 mapping = aes(x=variable, y=value, group = 1, shape=Repository))
  }

  if (min > min(silDf_melt$value)) {
    min = min(silDf_melt$value)
  }
  if (max < max(silDf_melt$value)) {
    max = max(silDf_melt$value)
  }
}

silPlot = silPlot +
  scale_x_continuous(name="k", breaks=1:5, labels=2:6) +
  scale_y_continuous(name="Quality", limits = c(min,max), breaks = c(0.25, 0.5, 0.7), labels=c(0.25, 0.5, 0.7)) +
  scale_colour_grey(start = 0.7, end = 0) +
  theme_calc()

#ggsave(plot = silPlot, filename=paste0(plotDir, "/silhouette_agro_obo_both.pdf"),
#       device="pdf", units="cm", width = 20, height = 10, dpi="print")

pg = plot_grid(stabPlot, silPlot, align = "v", nrow = 2, rel_heights = c(1/2, 1/2), labels=c("A", "B"))
#save_plot(paste0(plotDir, "/stability_silhouette_agro_obo_both.pdf"), pg, nrow=2, dpi="print")
#save_plot(paste0(plotDir, "/stability_silhouette_agro_obo_both.tiff"), pg, nrow=2, dpi="print", device="tiff")

#### CSV generation ----

# Stab
#write.csv(standardizeStabilityData(stabAgro, k.range), paste0(dataDir, "/stabilityAgro.csv"), row.names = TRUE)
#write.csv(standardizeStabilityData(stabObo, k.range), paste0(dataDir, "/stabilityObo.csv"), row.names = TRUE)
#write.csv(standardizeStabilityData(stabBoth, k.range), paste0(dataDir, "/stabilityAgro_Obo.csv"), row.names = TRUE)

# Qual
#write.csv(standardizeQualityData(qualAgro, k.range), paste0(dataDir, "/qualityAgro.csv"), row.names = TRUE)
#write.csv(standardizeQualityData(qualObo, k.range), paste0(dataDir, "/qualityObo.csv"), row.names = TRUE)
#write.csv(standardizeQualityData(qualBoth, k.range), paste0(dataDir, "/qualityAgro_Obo.csv"), row.names = TRUE)

#### Optimal k values ----
k.range=c(3,6)

# Agro
kOptAgro = getOptimalKValue(stabAgro, qualAgro, k.range=k.range)$Global_optimal_k
# Agro
kOptObo = getOptimalKValue(stabObo, qualObo, k.range=k.range)$Global_optimal_k
# Agro
kOptBoth = getOptimalKValue(stabBoth, qualBoth, k.range=k.range)$Global_optimal_k

tableKOpt = setNames(as.data.frame(matrix(nrow=19,ncol=3)), c("AgroPortal","OBO Foundry","AgroPortal + OBO Foundry"))
rownames(tableKOpt) <- meanStabAgro$Metric
tableKOpt["AgroPortal"] <- kOptAgro
tableKOpt["OBO Foundry"] <- kOptObo
tableKOpt["AgroPortal + OBO Foundry"] <- kOptBoth
tableKOpt

