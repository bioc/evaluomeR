library(evaluomeR)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(stringr)
library(reshape2)
library(ggthemes)

# Version
# Adding merged .tiff plots
# Neobernad committed on 24 Jan 2020

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
seed = 13000
plotDir=paste0(wd,"plots-rm")
dataDir=paste0(wd,"results-csv-rm")

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
csai16Path=paste0(wd,"data/data-csai16.csv")
csai17Path=paste0(wd,"data/data-csai17.csv")
csai18Path=paste0(wd,"data/data-csai18.csv")
csai19Path=paste0(wd,"data/data-csai19.csv")
csai20Path=paste0(wd,"data/data-csai20.csv")

csis16Path=paste0(wd,"data/data-csis16.csv")
csis17Path=paste0(wd,"data/data-csis17.csv")
csis18Path=paste0(wd,"data/data-csis18.csv")
csis19Path=paste0(wd,"data/data-csis19.csv")
csis20Path=paste0(wd,"data/data-csis20.csv")

orms16Path=paste0(wd,"data/data-orms16.csv")
orms17Path=paste0(wd,"data/data-orms17.csv")
orms18Path=paste0(wd,"data/data-orms18.csv")
orms19Path=paste0(wd,"data/data-orms19.csv")
orms20Path=paste0(wd,"data/data-orms20.csv")

inputCsai16 = read.csv(csai16Path, header = TRUE, sep = ";")
inputCsai17 = read.csv(csai17Path, header = TRUE, sep = ";")
inputCsai18 = read.csv(csai18Path, header = TRUE, sep = ";")
inputCsai19 = read.csv(csai19Path, header = TRUE, sep = ";")
inputCsai20 = read.csv(csai20Path, header = TRUE, sep = ";")

inputCsis16 = read.csv(csis16Path, header = TRUE, sep = ";")
inputCsis17 = read.csv(csis17Path, header = TRUE, sep = ";")
inputCsis18 = read.csv(csis18Path, header = TRUE, sep = ";")
inputCsis19 = read.csv(csis19Path, header = TRUE, sep = ";")
inputCsis20 = read.csv(csis20Path, header = TRUE, sep = ";")

inputOrms16 = read.csv(orms16Path, header = TRUE, sep = ";")
inputOrms17 = read.csv(orms17Path, header = TRUE, sep = ";")
inputOrms18 = read.csv(orms18Path, header = TRUE, sep = ";")
inputOrms19 = read.csv(orms19Path, header = TRUE, sep = ";")
inputOrms20 = read.csv(orms20Path, header = TRUE, sep = ";")

#### Stability [2,16] ----
k.range=c(2,15)
stabCsai16 <- stabilityRange(data=inputCsai16, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabCsai17 <- stabilityRange(data=inputCsai17, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabCsai18 <- stabilityRange(data=inputCsai18, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabCsai19 <- stabilityRange(data=inputCsai19, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabCsai20 <- stabilityRange(data=inputCsai20, k.range=k.range, bs=100, getImages = FALSE, seed=seed)

stabCsis16 <- stabilityRange(data=inputCsis16, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabCsis17 <- stabilityRange(data=inputCsis17, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabCsis18 <- stabilityRange(data=inputCsis18, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabCsis19 <- stabilityRange(data=inputCsis19, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabCsis20 <- stabilityRange(data=inputCsis20, k.range=k.range, bs=100, getImages = FALSE, seed=seed)

stabOrms16 <- stabilityRange(data=inputOrms16, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabOrms17 <- stabilityRange(data=inputOrms17, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabOrms18 <- stabilityRange(data=inputOrms18, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabOrms19 <- stabilityRange(data=inputOrms19, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabOrms20 <- stabilityRange(data=inputOrms20, k.range=k.range, bs=100, getImages = FALSE, seed=seed)

#### Stability SE to DF ----
meanStabCsai16 = standardizeStabilityData(stabCsai16, k.range)
meanStabCsai16$Metric = "CSAI16"
meanStabCsai16$category = "CSAI"
meanStabCsai16$year = "2016"
colnames(meanStabCsai16) = str_remove_all(colnames(meanStabCsai16), "k_")

meanStabCsai17 = standardizeStabilityData(stabCsai17, k.range)
meanStabCsai17$Metric = "CSAI17"
meanStabCsai17$category = "CSAI"
meanStabCsai17$year = "2017"
colnames(meanStabCsai17) = str_remove_all(colnames(meanStabCsai17), "k_")

meanStabCsai18 = standardizeStabilityData(stabCsai18, k.range)
meanStabCsai18$Metric = "CSAI18"
meanStabCsai18$category = "CSAI"
meanStabCsai18$year = "2018"
colnames(meanStabCsai18) = str_remove_all(colnames(meanStabCsai18), "k_")

meanStabCsai19 = standardizeStabilityData(stabCsai19, k.range)
meanStabCsai19$Metric = "CSAI19"
meanStabCsai19$category = "CSAI"
meanStabCsai19$year = "2019"
colnames(meanStabCsai19) = str_remove_all(colnames(meanStabCsai19), "k_")

meanStabCsai20 = standardizeStabilityData(stabCsai20, k.range)
meanStabCsai20$Metric = "CSAI20"
meanStabCsai20$category = "CSAI"
meanStabCsai20$year = "2020"
colnames(meanStabCsai20) = str_remove_all(colnames(meanStabCsai20), "k_")

meanStabCsis16 = standardizeStabilityData(stabCsis16, k.range)
meanStabCsis16$Metric = "CSIS16"
meanStabCsis16$category = "CSIS"
meanStabCsis16$year = "2016"
colnames(meanStabCsis16) = str_remove_all(colnames(meanStabCsis16), "k_")

meanStabCsis17 = standardizeStabilityData(stabCsis17, k.range)
meanStabCsis17$Metric = "CSIS17"
meanStabCsis17$category = "CSIS"
meanStabCsis17$year = "2017"
colnames(meanStabCsis17) = str_remove_all(colnames(meanStabCsis17), "k_")

meanStabCsis18 = standardizeStabilityData(stabCsis18, k.range)
meanStabCsis18$Metric = "CSIS18"
meanStabCsis18$category = "CSIS"
meanStabCsis18$year = "2018"
colnames(meanStabCsis18) = str_remove_all(colnames(meanStabCsis18), "k_")

meanStabCsis19 = standardizeStabilityData(stabCsis19, k.range)
meanStabCsis19$Metric = "CSIS19"
meanStabCsis19$category = "CSIS"
meanStabCsis19$year = "2019"
colnames(meanStabCsis19) = str_remove_all(colnames(meanStabCsis19), "k_")

meanStabCsis20 = standardizeStabilityData(stabCsis20, k.range)
meanStabCsis20$Metric = "CSIS20"
meanStabCsis20$category = "CSIS"
meanStabCsis20$year = "2020"
colnames(meanStabCsis20) = str_remove_all(colnames(meanStabCsis20), "k_")

meanStabOrms16 = standardizeStabilityData(stabOrms16, k.range)
meanStabOrms16$Metric = "ORMS16"
meanStabOrms16$category = "ORMS"
meanStabOrms16$year = "2016"
colnames(meanStabOrms16) = str_remove_all(colnames(meanStabOrms16), "k_")

meanStabOrms17 = standardizeStabilityData(stabOrms17, k.range)
meanStabOrms17$Metric = "ORMS17"
meanStabOrms17$category = "ORMS"
meanStabOrms17$year = "2017"
colnames(meanStabOrms17) = str_remove_all(colnames(meanStabOrms17), "k_")

meanStabOrms18 = standardizeStabilityData(stabOrms18, k.range)
meanStabOrms18$Metric = "ORMS18"
meanStabOrms18$category = "ORMS"
meanStabOrms18$year = "2018"
colnames(meanStabOrms18) = str_remove_all(colnames(meanStabOrms18), "k_")

meanStabOrms19 = standardizeStabilityData(stabOrms19, k.range)
meanStabOrms19$Metric = "ORMS19"
meanStabOrms19$category = "ORMS"
meanStabOrms19$year = "2019"
colnames(meanStabOrms19) = str_remove_all(colnames(meanStabOrms19), "k_")

meanStabOrms20 = standardizeStabilityData(stabOrms20, k.range)
meanStabOrms20$Metric = "ORMS20"
meanStabOrms20$category = "ORMS"
meanStabOrms20$year = "2020"
colnames(meanStabOrms20) = str_remove_all(colnames(meanStabOrms20), "k_")

#### Stability plotting (X k value, Y stability ----

stabDfList = list(meanStabCsai16, meanStabCsai17, meanStabCsai18, meanStabCsai19, meanStabCsai20,
                  meanStabCsis16, meanStabCsis17, meanStabCsis18, meanStabCsis19, meanStabCsis20,
                  meanStabOrms16, meanStabOrms17, meanStabOrms18, meanStabOrms19, meanStabOrms20)

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


#### Stability plotting by year (X k year, Y stability ----

max_stab_k_df = setNames(data.frame(matrix(ncol = 6, nrow = 0)),
                         c("Journal" ,"category", "year", "variable", "value", "ks"))
min = Inf
max = -Inf
for (stabDf in stabDfList) {
  stabDf[, "2"] = NULL
  stabDf_melt = melt(stabDf, id.vars = c("Metric", "category", "year"))
  colnames(stabDf_melt)[1] = "Journal"
  k_s = max(stabDf_melt$value)
  stabDf_melt = stabDf_melt[stabDf_melt$value==k_s, ]
  stabDf_melt$year = as.integer(stabDf_melt$year)

  # Knowing which k offers the highest stability (ks)
  for (colname in names(stabDf)) {
    if (! colname %in% c("Metric", "category", "year")) {
      if (stabDf[1, colname] == k_s) {
        stabDf_melt$ks = as.character(colname)
      }
    }
  }

  max_stab_k_df[nrow(max_stab_k_df) + 1, ] = stabDf_melt

  if (min > min(stabDf_melt$value)) {
    min = min(stabDf_melt$value)
  }
  if (max < max(stabDf_melt$value)) {
    max = max(stabDf_melt$value)
  }
}
max_stab_k_df

stabPlot_year = ggplot(data = max_stab_k_df) +
  #geom_text(mapping = aes(x=year, y=value,  label=k),  nudge_y = -0.01) +
  geom_point(mapping = aes(x=year, y=value,  shape=ks)) +
  geom_line(mapping = aes(x=year, y=value, linetype=category, colour=category, group=category))


stabPlot_year = stabPlot_year +
  scale_y_continuous(name="Stability", limits = c(0.75,1), breaks = c(0.6, 0.75, 0.85), labels=c(0.6, 0.75, 0.85)) +
  scale_colour_grey(start = 0.7, end = 0) +
  theme_calc(base_family = "sans")

ggsave(plot = stabPlot_year, filename=paste0(plotDir, "/stability_impact_factor_by_year.pdf"),
       device="pdf", units="cm", width = 20, height = 10, dpi="print")


#### Quality [3,15] ----
qualCsai16 <- qualityRange(data=inputCsai16, k.range=k.range, getImages = FALSE, seed=seed)
qualCsai17 <- qualityRange(data=inputCsai17, k.range=k.range, getImages = FALSE, seed=seed)
qualCsai18 <- qualityRange(data=inputCsai18, k.range=k.range, getImages = FALSE, seed=seed)
qualCsai19 <- qualityRange(data=inputCsai19, k.range=k.range, getImages = FALSE, seed=seed)
qualCsai20 <- qualityRange(data=inputCsai20, k.range=k.range, getImages = FALSE, seed=seed)

qualCsis16 <- qualityRange(data=inputCsis16, k.range=k.range, getImages = FALSE, seed=seed)
qualCsis17 <- qualityRange(data=inputCsis17, k.range=k.range, getImages = FALSE, seed=seed)
qualCsis18 <- qualityRange(data=inputCsis18, k.range=k.range, getImages = FALSE, seed=seed)
qualCsis19 <- qualityRange(data=inputCsis19, k.range=k.range, getImages = FALSE, seed=seed)
qualCsis20 <- qualityRange(data=inputCsis20, k.range=k.range, getImages = FALSE, seed=seed)

qualOrms16 <- qualityRange(data=inputOrms16, k.range=k.range, getImages = FALSE, seed=seed)
qualOrms17 <- qualityRange(data=inputOrms17, k.range=k.range, getImages = FALSE, seed=seed)
qualOrms18 <- qualityRange(data=inputOrms18, k.range=k.range, getImages = FALSE, seed=seed)
qualOrms19 <- qualityRange(data=inputOrms19, k.range=k.range, getImages = FALSE, seed=seed)
qualOrms20 <- qualityRange(data=inputOrms20, k.range=k.range, getImages = FALSE, seed=seed)

#### Quality SE to DF ----
silCsai16 = standardizeQualityData(qualCsai16, k.range)
silCsai16$Metric = "CSAI16"
silCsai16$category = "CSAI"
silCsai16$year = "2016"
colnames(silCsai16) = str_remove_all(colnames(silCsai16), "k_")

silCsai17 = standardizeQualityData(qualCsai17, k.range)
silCsai17$Metric = "CSAI17"
silCsai17$category = "CSAI"
silCsai17$year = "2017"
colnames(silCsai17) = str_remove_all(colnames(silCsai17), "k_")

silCsai18 = standardizeQualityData(qualCsai18, k.range)
silCsai18$Metric = "CSAI18"
silCsai18$category = "CSAI"
silCsai18$year = "2018"
colnames(silCsai18) = str_remove_all(colnames(silCsai18), "k_")

silCsai19 = standardizeQualityData(qualCsai19, k.range)
silCsai19$Metric = "CSAI19"
silCsai19$category = "CSAI"
silCsai19$year = "2019"
colnames(silCsai19) = str_remove_all(colnames(silCsai19), "k_")

silCsai20 = standardizeQualityData(qualCsai20, k.range)
silCsai20$Metric = "CSAI20"
silCsai20$category = "CSAI"
silCsai20$year = "2020"
colnames(silCsai20) = str_remove_all(colnames(silCsai20), "k_")

silCsis16 = standardizeQualityData(qualCsis16, k.range)
silCsis16$Metric = "CSIS16"
silCsis16$category = "CSIS"
silCsis16$year = "2016"
colnames(silCsis16) = str_remove_all(colnames(silCsis16), "k_")

silCsis17 = standardizeQualityData(qualCsis17, k.range)
silCsis17$Metric = "CSIS17"
silCsis17$category = "CSIS"
silCsis17$year = "2017"
colnames(silCsis17) = str_remove_all(colnames(silCsis17), "k_")

silCsis18 = standardizeQualityData(qualCsis18, k.range)
silCsis18$Metric = "CSIS18"
silCsis18$category = "CSIS"
silCsis18$year = "2018"
colnames(silCsis18) = str_remove_all(colnames(silCsis18), "k_")

silCsis19 = standardizeQualityData(qualCsis19, k.range)
silCsis19$Metric = "CSIS19"
silCsis19$category = "CSIS"
silCsis19$year = "2019"
colnames(silCsis19) = str_remove_all(colnames(silCsis19), "k_")

silCsis20 = standardizeQualityData(qualCsis20, k.range)
silCsis20$Metric = "CSIS20"
silCsis20$category = "CSIS"
silCsis20$year = "2020"
colnames(silCsis20) = str_remove_all(colnames(silCsis20), "k_")

silOrms16 = standardizeQualityData(qualOrms16, k.range)
silOrms16$Metric = "ORMS16"
silOrms16$category = "ORMS"
silOrms16$year = "2016"
colnames(silOrms16) = str_remove_all(colnames(silOrms16), "k_")

silOrms17 = standardizeQualityData(qualOrms17, k.range)
silOrms17$Metric = "ORMS17"
silOrms17$category = "ORMS"
silOrms17$year = "2017"
colnames(silOrms17) = str_remove_all(colnames(silOrms17), "k_")

silOrms18 = standardizeQualityData(qualOrms18, k.range)
silOrms18$Metric = "ORMS18"
silOrms18$category = "ORMS"
silOrms18$year = "2018"
colnames(silOrms18) = str_remove_all(colnames(silOrms18), "k_")

silOrms19 = standardizeQualityData(qualOrms19, k.range)
silOrms19$Metric = "ORMS19"
silOrms19$category = "ORMS"
silOrms19$year = "2019"
colnames(silOrms19) = str_remove_all(colnames(silOrms19), "k_")

silOrms20 = standardizeQualityData(qualOrms20, k.range)
silOrms20$Metric = "ORMS20"
silOrms20$category = "ORMS"
silOrms20$year = "2020"
colnames(silOrms20) = str_remove_all(colnames(silOrms20), "k_")

#### Quality plotting (X k value, Y quality ----

silDfList = list(silCsai16, silCsai17, silCsai18, silCsai19, silCsai20,
                 silCsis16, silCsis17, silCsis18, silCsis19, silCsis20,
                 silOrms16, silOrms17, silOrms18, silOrms19, silOrms20)

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


#### Quality plotting by year (X k year, Y quality ----

max_qual_k_df = setNames(data.frame(matrix(ncol = 6, nrow = 0)),
                         c("Journal" ,"category", "year", "variable", "value", "kg"))
min = Inf
max = -Inf
for (silDf in silDfList) {
  silDf[, "2"] = NULL
  silDf_melt = melt(silDf, id.vars = c("Metric", "category", "year"))
  colnames(silDf_melt)[1] = "Journal"
  k_g = max(silDf_melt$value)
  silDf_melt = silDf_melt[silDf_melt$value==max(silDf_melt$value), ]
  silDf_melt$year = as.integer(silDf_melt$year)


  # Knowing which k offers the highest goodness (kg)
  for (colname in names(silDf)) {
    if (! colname %in% c("Metric", "category", "year")) {
      if (silDf[1, colname] == k_g) {
        silDf_melt$kg = as.character(colname)
      }
    }
  }

  max_qual_k_df[nrow(max_qual_k_df) + 1, ] = silDf_melt

  if (min > min(silDf_melt$value)) {
    min = min(silDf_melt$value)
  }
  if (max < max(silDf_melt$value)) {
    max = max(silDf_melt$value)
  }
}

silPlot_year = ggplot(data = max_qual_k_df) +
  #geom_text(mapping = aes(x=year, y=value,  label=k),  nudge_y = -0.01) +
  geom_point(mapping = aes(x=year, y=value,  shape=kg)) +
  geom_line(mapping = aes(x=year, y=value, linetype=category, colour=category, group=category))


# limits = c(min, max)
silPlot_year = silPlot_year +
  scale_y_continuous(name="Quality", limits = c(0.5,0.7), breaks = c(0.25, 0.5, 0.7), labels=c(0.25, 0.5, 0.7)) +
  scale_colour_grey(start = 0.7, end = 0) +
  theme_calc()

ggsave(plot = stabPlot_year, filename=paste0(plotDir, "/silhouette_impact_factor_by_year.pdf"),
       device="pdf", units="cm", width = 20, height = 10, dpi="print")


# Saving A - B plot for stabPlot and silPlot
pg = plot_grid(stabPlot, silPlot, align = "v", nrow = 2, rel_heights = c(1/2, 1/2), labels=c("A", "B"))
save_plot(paste0(plotDir, "/stability_silhouette_impact_factor.pdf"), pg, nrow=2, dpi="print")
# save_plot(paste0(plotDir, "/stability_silhouette_impact_factor.png"), pg, nrow=2, dpi="print", device="png")

# Saving A - B plot for stabPlot_year and silPlot_year
pg = plot_grid(stabPlot_year, silPlot_year, align = "v", nrow = 2, rel_heights = c(1/2, 1/2), labels=c("A", "B"))
save_plot(paste0(plotDir, "/stability_silhouette_impact_factor_by_year.pdf"), pg, nrow=2, dpi="print")
# save_plot(paste0(plotDir, "/stability_silhouette_impact_factor_by_year.tiff"), pg, nrow=2, dpi="print", device="tiff")

#### CSV generation ----

# CSAI Stab
write.csv(standardizeStabilityData(stabCsai16, k.range), paste0(dataDir, "/stabilityCsai16.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabCsai17, k.range), paste0(dataDir, "/stabilityCsai17.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabCsai18, k.range), paste0(dataDir, "/stabilityCsai18.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabCsai19, k.range), paste0(dataDir, "/stabilityCsai19.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabCsai20, k.range), paste0(dataDir, "/stabilityCsai20.csv"), row.names = TRUE)

# CSAI Qual
write.csv(standardizeQualityData(qualCsai16, k.range), paste0(dataDir, "/qualityCsai16.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualCsai17, k.range), paste0(dataDir, "/qualityCsai17.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualCsai18, k.range), paste0(dataDir, "/qualityCsai18.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualCsai19, k.range), paste0(dataDir, "/qualityCsai19.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualCsai20, k.range), paste0(dataDir, "/qualityCsai20.csv"), row.names = TRUE)

# CSIS Stab
write.csv(standardizeStabilityData(stabCsis16, k.range), paste0(dataDir, "/stabilityCsis16.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabCsis17, k.range), paste0(dataDir, "/stabilityCsis17.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabCsis18, k.range), paste0(dataDir, "/stabilityCsis18.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabCsis19, k.range), paste0(dataDir, "/stabilityCsis19.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabCsis20, k.range), paste0(dataDir, "/stabilityCsis20.csv"), row.names = TRUE)

# CSIS Qual
write.csv(standardizeQualityData(qualCsis16, k.range), paste0(dataDir, "/qualityCsis16.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualCsis17, k.range), paste0(dataDir, "/qualityCsis17.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualCsis18, k.range), paste0(dataDir, "/qualityCsis18.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualCsis19, k.range), paste0(dataDir, "/qualityCsis19.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualCsis20, k.range), paste0(dataDir, "/qualityCsis20.csv"), row.names = TRUE)

# ORMS Stab
write.csv(standardizeStabilityData(stabOrms16, k.range), paste0(dataDir, "/stabilityOrms16.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabOrms17, k.range), paste0(dataDir, "/stabilityOrms17.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabOrms18, k.range), paste0(dataDir, "/stabilityOrms18.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabOrms19, k.range), paste0(dataDir, "/stabilityOrms19.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabOrms20, k.range), paste0(dataDir, "/stabilityOrms20.csv"), row.names = TRUE)

# ORMS Qual
write.csv(standardizeQualityData(qualOrms16, k.range), paste0(dataDir, "/qualityOrms16.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualOrms17, k.range), paste0(dataDir, "/qualityOrms17.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualOrms18, k.range), paste0(dataDir, "/qualityOrms18.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualOrms19, k.range), paste0(dataDir, "/qualityOrms19.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualOrms20, k.range), paste0(dataDir, "/qualityOrms20.csv"), row.names = TRUE)

#### Optimal k values ----
k.range=c(3,15)

# Csai
kOptCsai16 = getOptimalKValue(stabCsai16, qualCsai16, k.range=k.range)$Global_optimal_k
kOptCsai17 = getOptimalKValue(stabCsai17, qualCsai17, k.range=k.range)$Global_optimal_k
kOptCsai18 = getOptimalKValue(stabCsai18, qualCsai18, k.range=k.range)$Global_optimal_k
kOptCsai19 = getOptimalKValue(stabCsai19, qualCsai19, k.range=k.range)$Global_optimal_k
kOptCsai20 = getOptimalKValue(stabCsai20, qualCsai20, k.range=k.range)$Global_optimal_k

# Csis
kOptCsis16 = getOptimalKValue(stabCsis16, qualCsis16, k.range=k.range)$Global_optimal_k
kOptCsis17 = getOptimalKValue(stabCsis17, qualCsis17, k.range=k.range)$Global_optimal_k
kOptCsis18 = getOptimalKValue(stabCsis18, qualCsis18, k.range=k.range)$Global_optimal_k
kOptCsis19 = getOptimalKValue(stabCsis19, qualCsis19, k.range=k.range)$Global_optimal_k
kOptCsis20 = getOptimalKValue(stabCsis20, qualCsis20, k.range=k.range)$Global_optimal_k

# Orms
kOptOrms16 = getOptimalKValue(stabOrms16, qualOrms16, k.range=k.range)$Global_optimal_k
kOptOrms17 = getOptimalKValue(stabOrms17, qualOrms17, k.range=k.range)$Global_optimal_k
kOptOrms18 = getOptimalKValue(stabOrms18, qualOrms18, k.range=k.range)$Global_optimal_k
kOptOrms19 = getOptimalKValue(stabOrms19, qualOrms19, k.range=k.range)$Global_optimal_k
kOptOrms20 = getOptimalKValue(stabOrms20, qualOrms20, k.range=k.range)$Global_optimal_k

tableKOpt = setNames(as.data.frame(matrix(nrow=3,ncol=5)), c("2016","2017","2018", "2019", "2020"))
rownames(tableKOpt) <- c("CSAI", "CSIS", "ORMS")
tableKOpt["2016"] <- c(kOptCsai16, kOptCsis16, kOptOrms16)
tableKOpt["2017"] <- c(kOptCsai17, kOptCsis17, kOptOrms17)
tableKOpt["2018"] <- c(kOptCsai18, kOptCsis18, kOptOrms18)
tableKOpt["2019"] <- c(kOptCsai19, kOptCsis19, kOptOrms19)
tableKOpt["2020"] <- c(kOptCsai20, kOptCsis20, kOptOrms20)
tableKOpt
