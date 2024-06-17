library(evaluomeR)
library(DataCombine)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

evaluomeRSupportedCBI()
optimalKs=NULL
stabValues=NULL
silValues=NULL

cbi_used = c("kmeans","clara","clara_pam")
for (cbi in cbi_used) {
  GLOBAL_CBI = cbi
  optimalKs[[GLOBAL_CBI]] = NULL
  stabValues[[GLOBAL_CBI]] = NULL
  silValues[[GLOBAL_CBI]] = NULL
  source(paste0(wd,"agro.R"))
  source(paste0(wd,"obo.R"))
  optimalKs[[GLOBAL_CBI]]$Agro = as.data.frame(as.numeric(kOptTableAgro$Global_optimal_k),
                                               row.names=rownames(kOptTableObo))
  optimalKs[[GLOBAL_CBI]]$Obo = as.data.frame(as.numeric(kOptTableObo$Global_optimal_k),
                                              row.names=rownames(kOptTableObo))
  colnames(optimalKs[[GLOBAL_CBI]]$Agro) = c("k_opt")
  colnames(optimalKs[[GLOBAL_CBI]]$Obo) = c("k_opt")
  stabValues[[GLOBAL_CBI]]$Agro = stabilityTableAgro
  stabValues[[GLOBAL_CBI]]$Obo = stabilityTableObo
  silValues[[GLOBAL_CBI]]$Agro = silhouetteTableAgro
  silValues[[GLOBAL_CBI]]$Obo = silhouetteTableObo
}

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
outputDir=paste0(wd,"results-cbi")

#optimalKs = as.data.frame(optimalKs)
#rownames(optimalKs) = rownames(kOptTableObo)

# Data frame with metrics whose k value matches

#####################################################################################
# Table with stab./qual. scores of optimal k values for each CBI used
#####################################################################################


col_names = append(c("repository", "metric"), cbi_used)

stabilityComparison = setNames(data.frame(matrix(ncol = length(col_names), nrow = 0)), col_names)
qualityComparison = setNames(data.frame(matrix(ncol = length(col_names), nrow = 0)), col_names)
agroKComparison = setNames(data.frame(matrix(ncol = length(col_names), nrow = 0)), col_names)
oboKComparison = setNames(data.frame(matrix(ncol = length(col_names), nrow = 0)), col_names)
index = 1

for (metric in rownames(kOptTableObo)) {
  cbi_values = NULL
  for (cbi in cbi_used) { # Get stab. qual. values for each CBI for this metric
    agro_m_kopt = optimalKs[[cbi]]$Agro[metric, "k_opt"]
    obo_m_kopt = optimalKs[[cbi]]$Obo[metric, "k_opt"]
    agro_colname = paste0("k_", agro_m_kopt)
    obo_colname = paste0("k_", obo_m_kopt)
    agro_k_stability = as.numeric(stabValues[[cbi]]$Agro[metric, agro_colname])
    obo_k_stability = as.numeric(stabValues[[cbi]]$Obo[metric, obo_colname])
    agro_k_silhouette = as.numeric(silValues[[cbi]]$Agro[metric, agro_colname])
    obo_k_silhouette = as.numeric(silValues[[cbi]]$Obo[metric, obo_colname])
    cbi_values[[cbi]]$agro_k_stability = agro_k_stability
    cbi_values[[cbi]]$obo_k_stability = obo_k_stability
    cbi_values[[cbi]]$agro_k_silhouette = agro_k_silhouette
    cbi_values[[cbi]]$obo_k_silhouette = obo_k_silhouette
    cbi_values[[cbi]]$agro_k = agro_m_kopt
    cbi_values[[cbi]]$obo_k = obo_m_kopt
  }

  # Agro Stability
  new_row = c("Agroportal", metric,
              cbi_values[["kmeans"]]$agro_k_stability,
              cbi_values[["clara"]]$agro_k_stability,
              cbi_values[["clara_pam"]]$agro_k_stability
              )
  stabilityComparison <- InsertRow(stabilityComparison, NewRow = new_row, RowNum = index)

  # Agro Silhouette
  new_row = c("Agroportal", metric,
              cbi_values[["kmeans"]]$agro_k_silhouette,
              cbi_values[["clara"]]$agro_k_silhouette,
              cbi_values[["clara_pam"]]$agro_k_silhouette
              )
  qualityComparison <- InsertRow(qualityComparison, NewRow = new_row, RowNum = index)

  # Agro K Opt
  new_row = c("Agroportal", metric,
              cbi_values[["kmeans"]]$agro_k,
              cbi_values[["clara"]]$agro_k,
              cbi_values[["clara_pam"]]$agro_k
              )
  agroKComparison <- InsertRow(agroKComparison, NewRow = new_row, RowNum = index)

  # Obo Stability
  new_row = c("OBO Foundry", metric,
              cbi_values[["kmeans"]]$obo_k_stability,
              cbi_values[["clara"]]$obo_k_stability,
              cbi_values[["clara_pam"]]$obo_k_stability
              )
  stabilityComparison <- InsertRow(stabilityComparison, NewRow = new_row, RowNum = index)

  # Obo Silhouette
  new_row = c("OBO Foundry", metric,
              cbi_values[["kmeans"]]$obo_k_silhouette,
              cbi_values[["clara"]]$obo_k_silhouette,
              cbi_values[["clara_pam"]]$obo_k_silhouette
              )
  qualityComparison <- InsertRow(qualityComparison, NewRow = new_row, RowNum = index)

  # Obo K Opt
  new_row = c("Agroportal", metric,
              cbi_values[["kmeans"]]$obo_k,
              cbi_values[["clara"]]$obo_k,
              cbi_values[["clara_pam"]]$obo_k
  )
  oboKComparison <- InsertRow(oboKComparison, NewRow = new_row, RowNum = index)

  index = index + 1

}

stabilityComparison = na.omit(stabilityComparison)
qualityComparison = na.omit(qualityComparison)
agroKComparison = na.omit(agroKComparison)
oboKComparison = na.omit(oboKComparison)

stabilityComparison = stabilityComparison[order(stabilityComparison$metric), ]
qualityComparison = qualityComparison[order(qualityComparison$metric), ]
agroKComparison = agroKComparison[order(agroKComparison$metric), ]
oboKComparison = oboKComparison[order(oboKComparison$metric), ]


# Column name change from clara_pam to pam
colnames(stabilityComparison)[which(names(stabilityComparison) == "clara_pam")] <- "pam"
colnames(qualityComparison)[which(names(qualityComparison) == "clara_pam")] <- "pam"
colnames(agroKComparison)[which(names(agroKComparison) == "clara_pam")] <- "pam"
colnames(oboKComparison)[which(names(oboKComparison) == "clara_pam")] <- "pam"
measure.vars = c("kmeans", "clara", "pam")


# Stability Agro

agroStabilityComparison = stabilityComparison[stabilityComparison$repository=="Agroportal", ]
agroStabilityComparison$repository <- NULL
agroStabilityComparison = melt(agroStabilityComparison, id="metric", measure.vars = measure.vars)
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
  labs(title="Agroportal", subtitle="Stability scores with different algorithms",
       x="Metrics")

# Stability Obo

oboStabilityComparison = stabilityComparison[stabilityComparison$repository=="OBO Foundry", ]
oboStabilityComparison$repository <- NULL
oboStabilityComparison = melt(oboStabilityComparison, id="metric", measure.vars = measure.vars)
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
  labs(title="OBO Foundry", subtitle="Stability scores with different algorithms",
       x="Metrics")


# Silhouette Agro

agroSilhouetteComparison = qualityComparison[qualityComparison$repository=="Agroportal", ]
agroSilhouetteComparison$repository <- NULL
agroSilhouetteComparison = melt(agroSilhouetteComparison, id="metric", measure.vars = measure.vars)
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
  labs(title="Agroportal", subtitle="Goodness scores with different algorithms",
       x="Metrics")

# Silhouette Obo

oboSilhouetteComparison = qualityComparison[qualityComparison$repository=="OBO Foundry", ]
oboSilhouetteComparison$repository <- NULL
oboSilhouetteComparison = melt(oboSilhouetteComparison, id="metric", measure.vars = measure.vars)
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
  labs(title="OBO Foundry", subtitle="Goodness scores with different algorithms",
       x="Metrics")

# Huge plot!
#grid.arrange(agro_sta, obo_sta, agro_sil, obo_sil, ncol=1)
grid.arrange(agro_sta, agro_sil, ncol=1)
grid.arrange(obo_sta, obo_sil, ncol=1)


####################################################
# Dual plot with optimal K values across repositories and CBIs
####################################################

agroKComparisonMelt = agroKComparison
agroKComparisonMelt$repository <- NULL
agroKComparisonMelt = melt(agroKComparisonMelt, id="metric", measure.vars = measure.vars)
agroKComparisonMelt$value = as.numeric(agroKComparisonMelt$value)

agro_k <- ggplot(data=agroKComparisonMelt, aes(x=metric, y=value, group = variable)) +
  geom_line(aes(color=variable, linetype=variable)) +
  geom_point(aes(color=variable, shape=variable)) +
  geom_point(aes(color=variable)) +
  scale_colour_grey(start = 0, end = 0.5) +
  scale_y_continuous(name="K", breaks = c(2:15)) +
  coord_cartesian(ylim = c(2, 15)) + # For hollow shapes
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 16),
        legend.title = element_blank(),
        text = element_text(size=26),
        panel.grid.minor = element_blank()
  ) +
  labs(title="Agroportal", subtitle="Optimal K values with different algorithms",
       x="Metrics")


oboKComparisonMelt = oboKComparison
oboKComparisonMelt$repository <- NULL
oboKComparisonMelt = melt(oboKComparisonMelt, id="metric", measure.vars = measure.vars)
oboKComparisonMelt$value = as.numeric(oboKComparisonMelt$value)

obo_k <- ggplot(data=oboKComparisonMelt, aes(x=metric, y=value, group = variable)) +
  geom_line(aes(color=variable, linetype=variable)) +
  geom_point(aes(color=variable, shape=variable)) +
  geom_point(aes(color=variable)) +
  scale_colour_grey(start = 0, end = 0.5) +
  scale_y_continuous(name="K", breaks = c(2:15)) +
  coord_cartesian(ylim = c(2, 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 16),
        legend.title = element_blank(),
        text = element_text(size=26),
        panel.grid.minor = element_blank()
  ) +
  labs(title="OBO Foundry", subtitle="Optimal K values with different algorithms",
       x="Metrics")

# Plots
agro_k
obo_k
grid.arrange(agro_k, obo_k, ncol=1)

####################################################
# Table with the classification of each stab./qual. score and CBI
####################################################

# Agro stab.
agroKOptStaLabels = stabilityComparison[stabilityComparison$repository=="Agroportal", ]
agroKOptStaLabelsCopy = agroKOptStaLabels

for (cbi in measure.vars) {
  agroKOptStaLabels[
    which(agroKOptStaLabelsCopy[cbi] > 0.85 & agroKOptStaLabelsCopy[cbi] <= 1), cbi
    ] = "Highly stable"
  agroKOptStaLabels[
    which(agroKOptStaLabelsCopy[cbi] > 0.75 & agroKOptStaLabelsCopy[cbi] <= 0.85), cbi
    ] = "Stable"
  agroKOptStaLabels[
    which(agroKOptStaLabelsCopy[cbi] >= 0.60 & agroKOptStaLabelsCopy[cbi] <= 0.75), cbi
    ] = "Doubtful"
  agroKOptStaLabels[
    which(agroKOptStaLabelsCopy[cbi] >= 0 & agroKOptStaLabelsCopy[cbi] < 0.60), cbi
    ] = "Unstable"
}

# Agro sil.
agroKOptSilLabels = qualityComparison[qualityComparison$repository=="Agroportal", ]
agroKOptSilLabelsCopy = agroKOptSilLabels

for (cbi in measure.vars) {
  agroKOptSilLabels[
    which(agroKOptSilLabelsCopy[cbi] > 0.7 & agroKOptSilLabelsCopy[cbi] <= 1), cbi
    ] = "Strong clust. struct."
  agroKOptSilLabels[
    which(agroKOptSilLabelsCopy[cbi] > 0.5 & agroKOptSilLabelsCopy[cbi] <= 0.7), cbi
    ] = "Reasonable clust. struct."
  agroKOptSilLabels[
    which(agroKOptSilLabelsCopy[cbi] > 0.25 & agroKOptSilLabelsCopy[cbi] <= 0.50), cbi
    ] = "Weak clust. struct."
  agroKOptSilLabels[
    which(agroKOptSilLabelsCopy[cbi] >= -1 & agroKOptSilLabelsCopy[cbi] <= 0.25), cbi
    ] = "No substantial clust. struct."
}

# OBO stab.
oboKOptStaLabels = stabilityComparison[stabilityComparison$repository=="OBO Foundry", ]
oboKOptStaLabelsCopy = oboKOptStaLabels

for (cbi in measure.vars) {
  oboKOptStaLabels[
    which(oboKOptStaLabelsCopy[cbi] > 0.85 & oboKOptStaLabelsCopy[cbi] <= 1), cbi
    ] = "Highly stable"
  oboKOptStaLabels[
    which(oboKOptStaLabelsCopy[cbi] > 0.75 & oboKOptStaLabelsCopy[cbi] <= 0.85), cbi
    ] = "Stable"
  oboKOptStaLabels[
    which(oboKOptStaLabelsCopy[cbi] >= 0.60 & oboKOptStaLabelsCopy[cbi] <= 0.75), cbi
    ] = "Doubtful"
  oboKOptStaLabels[
    which(oboKOptStaLabelsCopy[cbi] >= 0 & oboKOptStaLabelsCopy[cbi] < 0.60), cbi
    ] = "Unstable"
}

# OBO sil.
oboKOptSilLabels = qualityComparison[qualityComparison$repository=="OBO Foundry", ]
oboKOptSilLabelsCopy = oboKOptSilLabels

for (cbi in measure.vars) {
  oboKOptSilLabels[
    which(oboKOptSilLabelsCopy[cbi] > 0.7 & oboKOptSilLabelsCopy[cbi] <= 1), cbi
    ] = "Strong clust. struct."
  oboKOptSilLabels[
    which(oboKOptSilLabelsCopy[cbi] > 0.5 & oboKOptSilLabelsCopy[cbi] <= 0.7), cbi
    ] = "Reasonable clust. struct."
  oboKOptSilLabels[
    which(oboKOptSilLabelsCopy[cbi] > 0.25 & oboKOptSilLabelsCopy[cbi] <= 0.50), cbi
    ] = "Weak clust. struct."
  oboKOptSilLabels[
    which(oboKOptSilLabelsCopy[cbi] >= -1 & oboKOptSilLabelsCopy[cbi] <= 0.25), cbi
    ] = "No substantial clust. struct."
}

# Export

for (cbi in measure.vars) {
  kOptLabelsTable = NULL
  kOptLabelsTable$agro_stab = agroKOptStaLabels[cbi]
  kOptLabelsTable$agro_sil = agroKOptSilLabels[cbi]
  kOptLabelsTable$obo_stab = oboKOptStaLabels[cbi]
  kOptLabelsTable$obo_sil = oboKOptSilLabels[cbi]
  kOptLabelsTable = as.data.frame(kOptLabelsTable)
  rownames(kOptLabelsTable) = oboKOptSilLabels$metric
  colnames(kOptLabelsTable) = c(paste0("Agro_stab_", cbi), paste0("Agro_sil_", cbi),
                                paste0("Obo_stab_", cbi), paste0("Obo_sil_", cbi))
  csvPath = paste0(outputDir, "/kopt_table_cuanlitative_",cbi,".csv")
  write.csv(kOptLabelsTable, csvPath, row.names = TRUE)
}


