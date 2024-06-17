library(evaluomeR)

individuals_per_cluster = function(qualityResult) {
  qual_df = as.data.frame(assay(qualityResult))


  cluster_pos_str = as.character(unlist(qual_df["Cluster_position"]))
  cluster_labels_str = as.character(unlist(qual_df["Cluster_labels"]))

  cluster_pos = as.list(strsplit(cluster_pos_str, ",")[[1]])
  cluster_labels = as.list(strsplit(cluster_labels_str, ",")[[1]])

  individuals_in_cluster = as.data.frame(cbind(cluster_labels, cluster_pos))
  colnames(individuals_in_cluster) = c("Individual", "InCluster")

  return(individuals_in_cluster)
}

data("ontMetrics")
metricsRelevancy = getMetricsRelevancy(ontMetrics, k=3, alpha=0.1, seed=100)
# RSKC output object
metricsRelevancy$rskc
# Trimmed cases from input (row indexes)
metricsRelevancy$trimmed_cases
# Metrics relevancy table
metricsRelevancy$relevancy


test = qualityRange(data=ontMetrics, k.range=c(3,3),
                             seed=13007,
                             all_metrics=TRUE,
                             cbi="rskc", L1=2, alpha=0)

# Shows how clusters are partitioned according to the individuals
individuals_per_cluster(test$k_3)
