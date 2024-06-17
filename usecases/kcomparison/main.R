library(evaluomeR)
library(DataCombine)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")

evaluomeRSupportedCBI()
GLOBAL_CBI = "kmeans"

source(paste0(wd,"agro.R"))
source(paste0(wd,"obo.R"))
outputDir=paste0(wd,"results-k")

# Table with optimal k values for Agro and OBO
optimalKs=NULL
optimalKs$Agro = as.numeric(kOptTableAgro$Global_optimal_k)
optimalKs$Obo = as.numeric(kOptTableObo$Global_optimal_k)
optimalKs = as.data.frame(optimalKs)
rownames(optimalKs) = rownames(kOptTableObo)


# Table with suboptimal and optimal k values for Agro and OBO
table = NULL
table$Agro_stability_max_k = as.integer(kOptTableAgro$Stability_max_k)
table$Agro_quality_max_k = as.integer(kOptTableAgro$Quality_max_k)
table$Agro_global_max_k = as.integer(kOptTableAgro$Global_optimal_k)

table$Obo_stability_max_k = as.integer(kOptTableObo$Stability_max_k)
table$Obo_quality_max_k = as.integer(kOptTableObo$Quality_max_k)
table$Obo_global_max_k = as.integer(kOptTableObo$Global_optimal_k)

table = as.data.frame(table)
rownames(table) = rownames(kOptTableObo)

csvPath = paste0(outputDir, "/table_", GLOBAL_CBI,".csv")
#write.csv(table, csvPath, row.names = TRUE)

# Data frame with metrics whose k value matches
sameKIndexes = which(optimalKs$Agro == optimalKs$Obo)
sameKMetrics = rownames(optimalKs)[sameKIndexes]
sameKDf = optimalKs[sameKMetrics,]

#####################################################################################
# Table with stab./qual. scores of optimal k values vs K=5
#####################################################################################

stabilityComparison = setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("repository", "metric", "k opt.", "k=5"))
qualityComparison = setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("repository", "metric", "k opt.", "k=5"))
stabCount = 1
silCount = 1

for (metric in rownames(table)) {
  agro_m_kopt = table[metric, "Agro_global_max_k"]
  obo_m_kopt = table[metric, "Obo_global_max_k"]
  agro_colname = paste0("k_", agro_m_kopt)
  obo_colname = paste0("k_", obo_m_kopt)
  agro_k_stability = as.numeric(stabilityTableAgro[metric, agro_colname])
  obo_k_stability = as.numeric(stabilityTableObo[metric, obo_colname])
  agro_k_silhouette = as.numeric(silhouetteTableAgro[metric, agro_colname])
  obo_k_silhouette = as.numeric(silhouetteTableObo[metric, obo_colname])

  # Agro Stability
  new_row = c("Agroportal", metric, as.numeric(agro_k_stability), stabilityTableAgro[metric, "k_5"])
  stabilityComparison <- InsertRow(stabilityComparison, NewRow = new_row, RowNum = stabCount)
  stabCount = stabCount + 1

  # Agro Silhouette
  new_row = c("Agroportal", metric, agro_k_silhouette, silhouetteTableAgro[metric, "k_5"])
  qualityComparison <- InsertRow(qualityComparison, NewRow = new_row, RowNum = silCount)
  silCount = silCount + 1

  # Obo Stability
  new_row = c("OBO Foundry", metric, obo_k_stability, stabilityTableObo[metric, "k_5"])
  stabilityComparison <- InsertRow(stabilityComparison, NewRow = new_row, RowNum = stabCount)
  stabCount = stabCount + 1

  # Obo Silhouette
  new_row = c("OBO Foundry", metric, obo_k_silhouette, silhouetteTableObo[metric, "k_5"])
  qualityComparison <- InsertRow(qualityComparison, NewRow = new_row, RowNum = silCount)
  silCount = silCount + 1

}

stabilityComparison = na.omit(stabilityComparison)
qualityComparison = na.omit(qualityComparison)

kOptTable = NULL
kOptTable$agro_stab = stabilityComparison[stabilityComparison$repository == "Agroportal", ]$`k opt.`
kOptTable$agro_sil = qualityComparison[qualityComparison$repository == "Agroportal", ]$`k opt.`
kOptTable$obo_stab = stabilityComparison[stabilityComparison$repository == "OBO Foundry", ]$`k opt.`
kOptTable$obo_sil = qualityComparison[qualityComparison$repository == "OBO Foundry", ]$`k opt.`
kOptTable = as.data.frame(kOptTable)
rownames(kOptTable) = stabilityComparison[stabilityComparison$repository=="Agroportal", ]$metric

csvPath = paste0(outputDir, "/stability_qualityKopt",".csv")
#write.csv(kOptTable, csvPath, row.names = TRUE)

# Stability Agro

agroStabilityComparison = stabilityComparison[stabilityComparison$repository=="Agroportal", ]
agroStabilityComparison$repository <- NULL
agroStabilityComparison = melt(agroStabilityComparison, id="metric", measure.vars = c("k opt.", "k=5"))
agroStabilityComparison$value = as.numeric(agroStabilityComparison$value)

agro_sta <- ggplot(data=agroStabilityComparison, aes(x=metric, y=value, group = variable)) +
  geom_line(aes(color=variable, linetype=variable)) +
  geom_point(aes(color=variable, shape=variable)) +
  geom_point(aes(color=variable)) +
  scale_colour_grey(start = 0, end = 0.5) +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_y_continuous(name="Stability") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.title = element_blank(),
        text = element_text(size=26)
  ) +
  labs(title="Agroportal", subtitle="Stability scores for k opt. and k=5",
       x="Metrics")

# Stability Obo

oboStabilityComparison = stabilityComparison[stabilityComparison$repository=="OBO Foundry", ]
oboStabilityComparison$repository <- NULL
oboStabilityComparison = melt(oboStabilityComparison, id="metric", measure.vars = c("k opt.", "k=5"))
oboStabilityComparison$value = as.numeric(oboStabilityComparison$value)

obo_sta <- ggplot(data=oboStabilityComparison, aes(x=metric, y=value, group = variable)) +
  geom_line(aes(color=variable, linetype=variable)) +
  geom_point(aes(color=variable, shape=variable)) +
  geom_point(aes(color=variable)) +
  scale_colour_grey(start = 0, end = 0.5) +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_y_continuous(name="Stability") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.title = element_blank(),
        text = element_text(size=26)
  ) +
  labs(title="OBO Foundry", subtitle="Stability scores for k opt. and k=5",
       x="Metrics")


# Silhouette Agro

agroSilhouetteComparison = qualityComparison[qualityComparison$repository=="Agroportal", ]
agroSilhouetteComparison$repository <- NULL
agroSilhouetteComparison = melt(agroSilhouetteComparison, id="metric", measure.vars = c("k opt.", "k=5"))
agroSilhouetteComparison$value = as.numeric(agroSilhouetteComparison$value)

agro_sil <- ggplot(data=agroSilhouetteComparison, aes(x=metric, y=value, group = variable)) +
  geom_line(aes(color=variable, linetype=variable)) +
  geom_point(aes(color=variable, shape=variable)) +
  geom_point(aes(color=variable)) +
  scale_colour_grey(start = 0, end = 0.5) +
  scale_y_continuous(name="Goodness") +
  coord_cartesian(ylim = c(0.5, 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.title = element_blank(),
        text = element_text(size=26)
  ) +
  labs(title="Agroportal", subtitle="Goodness scores for k opt. and k=5",
       x="Metrics")

# Silhouette Obo

oboSilhouetteComparison = qualityComparison[qualityComparison$repository=="OBO Foundry", ]
oboSilhouetteComparison$repository <- NULL
oboSilhouetteComparison = melt(oboSilhouetteComparison, id="metric", measure.vars = c("k opt.", "k=5"))
oboSilhouetteComparison$value = as.numeric(oboSilhouetteComparison$value)

obo_sil <- ggplot(data=oboSilhouetteComparison, aes(x=metric, y=value, group = variable)) +
  geom_line(aes(color=variable, linetype=variable)) +
  geom_point(aes(color=variable, shape=variable)) +
  geom_point(aes(color=variable)) +
  scale_colour_grey(start = 0, end = 0.5) +
  scale_y_continuous(name="Goodness") + # col="# of cylinders"
  coord_cartesian(ylim = c(0.5, 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.title = element_blank(),
        text = element_text(size=26)
  ) +
  labs(title="OBO Foundry", subtitle="Goodness scores for k opt. and k=5",
       x="Metrics")

# Huge plot!
grid.arrange(agro_sta, obo_sta, agro_sil, obo_sil, ncol=1)

#####################################################################################
# Table with stab./qual. scores of optimal k in AgroPortal and OBo Foundry
#####################################################################################

# Stability Agro

agroStabilityComparison = stabilityComparison[stabilityComparison$repository=="Agroportal", ]
agroStabilityComparison$`k=5` <- NULL
agroStabilityComparison = melt(agroStabilityComparison, id="metric", measure.vars = c("k opt."))
agroStabilityComparison$value = as.numeric(agroStabilityComparison$value)
agroStabilityComparison$variable = rep("AgroPortal", length(agroStabilityComparison$variable))

# Stability OBO

oboStabilityComparison = stabilityComparison[stabilityComparison$repository=="OBO Foundry", ]
oboStabilityComparison$`k=5` <- NULL
oboStabilityComparison = melt(oboStabilityComparison, id="metric", measure.vars = c("k opt."))
oboStabilityComparison$value = as.numeric(oboStabilityComparison$value)
oboStabilityComparison$variable = rep("OBO Foundry", length(oboStabilityComparison$variable))

toPlot = rbind(agroStabilityComparison, oboStabilityComparison)

stabPlot = ggplot(data=toPlot, aes(x=metric, y=value, group = variable)) +
  geom_line(aes(color=variable, linetype=variable)) +
  geom_point(aes(color=variable, shape=variable)) +
  scale_colour_grey(start = 0, end = 0.5) +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_y_continuous(name="Stability") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.title = element_blank(),
        text = element_text(size=26)
  ) +
  labs(title=NULL, subtitle=NULL,
       x="Metrics")

# Goodness Agro

agroSilhouetteComparison = qualityComparison[qualityComparison$repository=="Agroportal", ]
agroSilhouetteComparison$`k=5` <- NULL
agroSilhouetteComparison = melt(agroSilhouetteComparison, id="metric", measure.vars = c("k opt."))
agroSilhouetteComparison$value = as.numeric(agroSilhouetteComparison$value)
agroSilhouetteComparison$variable = rep("AgroPortal", length(agroStabilityComparison$variable))

# Goodness OBO

oboSilhouetteComparison = qualityComparison[qualityComparison$repository=="OBO Foundry", ]
oboSilhouetteComparison$`k=5` <- NULL
oboSilhouetteComparison = melt(oboSilhouetteComparison, id="metric", measure.vars = c("k opt."))
oboSilhouetteComparison$value = as.numeric(oboSilhouetteComparison$value)
oboSilhouetteComparison$variable = rep("OBO Foundry", length(oboStabilityComparison$variable))

toPlot = rbind(agroSilhouetteComparison, oboSilhouetteComparison)

goodnessPlot <- ggplot(data=toPlot, aes(x=metric, y=value, group = variable)) +
  geom_line(aes(color=variable, linetype=variable)) +
  geom_point(aes(color=variable, shape=variable)) +
  geom_point(aes(color=variable)) +
  scale_colour_grey(start = 0, end = 0.5) +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_y_continuous(name="Goodness") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.title = element_blank(),
        text = element_text(size=26)
  ) +
  labs(title=NULL, subtitle=NULL,
       x="Metrics")

grid.arrange(stabPlot, goodnessPlot, ncol=1)

####################################################
# Table with the classification of each stab./qual. score
####################################################

optimal_k = "k opt."
kOptLabelsTable = NULL

# Agro stab.
agroKOptStaLabels = stabilityComparison[stabilityComparison$repository=="Agroportal", ]
agroKOptStaLabels[, c("repository", "k=5")] <- list(NULL)
agroKOptStaLabelsCopy = agroKOptStaLabels

agroKOptStaLabels[
  which(agroKOptStaLabelsCopy[optimal_k] > 0.85 & agroKOptStaLabelsCopy[optimal_k] <= 1), optimal_k
  ] = "Highly stable"
agroKOptStaLabels[
  which(agroKOptStaLabelsCopy[optimal_k] > 0.75 & agroKOptStaLabelsCopy[optimal_k] <= 0.85), optimal_k
  ] = "Stable"
agroKOptStaLabels[
  which(agroKOptStaLabelsCopy[optimal_k] >= 0.60 & agroKOptStaLabelsCopy[optimal_k] <= 0.75), optimal_k
  ] = "Doubtful"
agroKOptStaLabels[
  which(agroKOptStaLabelsCopy[optimal_k] >= 0 & agroKOptStaLabelsCopy[optimal_k] < 0.60), optimal_k
  ] = "Unstable"

# Agro sil.
agroKOptSilLabels = qualityComparison[qualityComparison$repository=="Agroportal", ]
agroKOptSilLabels[, c("repository", "k=5")] <- list(NULL)
agroKOptSilLabelsCopy = agroKOptSilLabels

agroKOptSilLabels[
  which(agroKOptSilLabelsCopy[optimal_k] > 0.7 & agroKOptSilLabelsCopy[optimal_k] <= 1), optimal_k
  ] = "Strong clust. struct."
agroKOptSilLabels[
  which(agroKOptSilLabelsCopy[optimal_k] > 0.5 & agroKOptSilLabelsCopy[optimal_k] <= 0.7), optimal_k
  ] = "Reasonable clust. struct."
agroKOptSilLabels[
  which(agroKOptSilLabelsCopy[optimal_k] > 0.25 & agroKOptSilLabelsCopy[optimal_k] <= 0.50), optimal_k
  ] = "Weak clust. struct."
agroKOptSilLabels[
  which(agroKOptSilLabelsCopy[optimal_k] >= -1 & agroKOptSilLabelsCopy[optimal_k] <= 0.25), optimal_k
  ] = "No substantial clust. struct."

# OBO stab.
oboKOptStaLabels = stabilityComparison[stabilityComparison$repository=="OBO Foundry", ]
oboKOptStaLabels[, c("repository", "k=5")] <- list(NULL)
oboKOptStaLabelsCopy = oboKOptStaLabels

oboKOptStaLabels[
  which(oboKOptStaLabelsCopy[optimal_k] > 0.85 & oboKOptStaLabelsCopy[optimal_k] <= 1), optimal_k
  ] = "Highly stable"
oboKOptStaLabels[
  which(oboKOptStaLabelsCopy[optimal_k] > 0.75 & oboKOptStaLabelsCopy[optimal_k] <= 0.85), optimal_k
  ] = "Stable"
oboKOptStaLabels[
  which(oboKOptStaLabelsCopy[optimal_k] >= 0.60 & oboKOptStaLabelsCopy[optimal_k] <= 0.75), optimal_k
  ] = "Doubtful"
oboKOptStaLabels[
  which(oboKOptStaLabelsCopy[optimal_k] >= 0 & oboKOptStaLabelsCopy[optimal_k] < 0.60), optimal_k
  ] = "Unstable"

# OBO sil.
oboKOptSilLabels = qualityComparison[qualityComparison$repository=="OBO Foundry", ]
oboKOptSilLabels[, c("repository", "k=5")] <- list(NULL)
oboKOptSilLabelsCopy = oboKOptSilLabels

oboKOptSilLabels[
  which(oboKOptSilLabelsCopy[optimal_k] > 0.7 & oboKOptSilLabelsCopy[optimal_k] <= 1), optimal_k
  ] = "Strong clust. struct."
oboKOptSilLabels[
  which(oboKOptSilLabelsCopy[optimal_k] > 0.5 & oboKOptSilLabelsCopy[optimal_k] <= 0.7), optimal_k
  ] = "Reasonable clust. struct."
oboKOptSilLabels[
  which(oboKOptSilLabelsCopy[optimal_k] > 0.25 & oboKOptSilLabelsCopy[optimal_k] <= 0.50), optimal_k
  ] = "Weak clust. struct."
oboKOptSilLabels[
  which(oboKOptSilLabelsCopy[optimal_k] >= -1 & oboKOptSilLabelsCopy[optimal_k] <= 0.25), optimal_k
  ] = "No substantial clust. struct."


kOptLabelsTable$agro_stab = agroKOptStaLabels$`k opt.`
kOptLabelsTable$agro_sil = agroKOptSilLabels$`k opt.`
kOptLabelsTable$obo_stab = oboKOptStaLabels$`k opt.`
kOptLabelsTable$obo_sil = oboKOptSilLabels$`k opt.`
kOptLabelsTable = as.data.frame(kOptLabelsTable)
rownames(kOptLabelsTable) = stabilityComparison[stabilityComparison$repository=="Agroportal", ]$metric
csvPath = paste0(outputDir, "/both_repo_k_values_cualitative",".csv")
#write.csv(kOptLabelsTable, csvPath, row.names = TRUE)

length(which(stabilityComparison[stabilityComparison$repository=="Agroportal", "k opt."] > 0.75))
length(which(stabilityComparison[stabilityComparison$repository=="OBO Foundry", "k opt."] > 0.75))
length(which(qualityComparison[qualityComparison$repository=="Agroportal", "k opt."] > 0.5))
length(which(qualityComparison[qualityComparison$repository=="OBO Foundry", "k opt."] > 0.5))

##
# Plotting signatures of Agro and OBO for their optimal k value vs k = 5
##

agroVector = as.integer(kOptTableAgro$Global_optimal_k)
oboVector = as.integer(kOptTableObo$Global_optimal_k)
k5Vector = rep(5, length(agroVector))
plotMetricsClusterComparison(inputDataAgro, agroVector, k5Vector, seed = 13606)
plotMetricsClusterComparison(inputDataObo, oboVector, k5Vector, seed = 13606)
