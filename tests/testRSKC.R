library(evaluomeR)
library(RSKC)
library(sparcl)


# Dataframe for the use case is 'ontMetrics' provided by our evaluomeR package.

data("ontMetrics")
df = as.data.frame(assay(ontMetrics))
df["Description"] = NULL # Description column not relevant atm.
head(df, 5)
data("ontMetrics")

# RSKC
# Robust and Sparse K-Means clustering [[1]](#1) requires to select mainly three parameters:
# - **nlc**: Number of *K* cluster. It is to be determined by *evaluomeR* optimal *K* algorithm.
#   - **L<sub>1</sub>**: The tuning parameter for sparce clustering. It acts as the upper bound restraint for the vector of weights. 1 $<$ L<sub>1</sub> $\leq$ $\sqrt{num.variables}$.
#   - **$\alpha$**: The trimming portion [[4]](#4) used in the robust clustering.

# Optimal K clusters value
# Here, we make use of *evaluomeR* to figure out the optimal $k$ value. The algorithm on how the optimal is calculated is outlined in [[7]](#7). We consider the $k$ range [3,15] for the analysis of the optimal $k$, avoiding $k=2$ to prevent from having binary classifications.
seed=100
k.range=c(3,15)
stabilityData <- stabilityRange(data=ontMetrics, k.range=k.range, bs=20, getImages = FALSE, seed=seed)
qualityData <- qualityRange(data=ontMetrics, k.range=k.range, getImages = FALSE, seed=seed)
optK <- getOptimalKValue(stabilityData, qualityData, k.range=k.range)

# Optimal $k$ values individually per input metric are:
optK[c('Metric','Global_optimal_k')]

k_values = as.numeric(unlist(optK['Global_optimal_k']))
global_k_value = floor(mean(k_values))
print(paste0("Taking global optimal K value: ", global_k_value))

plotMetricsClusterComparison(ontMetrics, k.vector1=global_k_value)

# Figuring out the L1 upper boundry
# In [[2]](#2) authors provide description of the algorithm to select the tunning parameter L<sub>1</sub> for
# the sparse K-means, which consist of independent permutations from the same source data matrix and the gap
# statistic [[5]](#5). This algorithm for tuning the L<sub>1</sub> parameter and others described in [[2]](#2)
# are presented in 'sparcl' R package [[6]](#6).

dataMatrix = as.matrix(df)
dataMatrix = scale(dataMatrix, TRUE, TRUE)
head(dataMatrix, 5)

# Considering that for the dataset the global optimal $k$ is $k=4$, we can now compute the
# permutations to figure out the boundry L<sub>1</sub> with the method 'KMeansSparseCluster.permute'
# from 'sparcl' [[6]](#6).

# Note: 1 $<$ L<sub>1</sub> $\leq$ $\sqrt{num.variables}$.

wbounds = seq(2,sqrt(ncol(dataMatrix)), len=30)
km.perm <- KMeansSparseCluster.permute(dataMatrix,K=global_k_value,wbounds=wbounds,nperms=5)
print(km.perm)
plot(km.perm)

l1 = km.perm$bestw
print(paste0("Best L1 upper bound is: ", l1))



# Metrics relevancy
rskc_out = RSKC(df["ANOnto"], global_k_value, 0.1, L1 = l1, nstart = 200,
                silent=TRUE, scaling = FALSE, correlation = FALSE)
cat(paste0("L1 value: ", l1,"\n"))
cat(names(rskc_out$weights)[1], ": ", rskc_out$weights[1],"\n")
cat(names(rskc_out$weights)[2], ": ", rskc_out$weights[2],"\n")
cat(names(rskc_out$weights)[3], ": ", rskc_out$weights[3],"\n")
cat("---\n")

rskc_out

# Trimmed cases:

# oE: Indices of the cases trimmed in squared Euclidean distances.
# oW:	Indices of the cases trimmed in weighted squared Euclidean distances. If L1 =NULL,
# then oW are the cases trimmed in the Euclidean distance, because all the features have the same weights, i.e., 1's.
union_vector = c(rskc_out$oE,rskc_out$oW)
union_vector_unique = unique(union_vector)
union_vector_unique = sort(union_vector_unique)

print(paste0("Trimmed cases from input dataframe: "))
union_vector_unique

options(scipen=10)

columns = c('metric', 'weight')
rskc_df = data.frame(matrix(ncol = length(columns), nrow = length(rskc_out$weights)))
colnames(rskc_df) = columns
rskc_df['metric'] = names(rskc_out$weights)
rskc_df['weight'] = rskc_out$weights
rskc_df

# Relevancy table
rskc_df_sorted = rskc_df[order(rskc_df$weight, decreasing = TRUE), ]
rskc_df_sorted

# References <a class="anchor" id="references"></a>

#<a id="1">[1]</a>
#  Kondo, Y., Salibian-Barrera, M., & Zamar, R. (2016). RSKC: An R Package for a Robust and Sparse K-Means Clustering Algorithm. Journal of Statistical Software, 72(5), 1–26. https://doi.org/10.18637/jss.v072.i05

#<a id="2">[2]</a>
#  Witten, D. M., & Tibshirani, R. (2010). A framework for feature selection in clustering. Journal of the American Statistical Association, 105(490), 713–726. https://doi.org/10.1198/jasa.2010.tm09415

#<a id="3">[3]</a>
#  Robert Tibshirani, & Guenther Walther (2005). Cluster Validation by Prediction Strength. Journal of Computational and Graphical Statistics, 14(3), 511-528. https://doi.org/10.1198/106186005X59243

#<a id="4">[4]</a>
#  Gordaliza, A. (1991). On the breakdown point of multivariate location estimators based on trimming procedures. Statistics & Probability Letters, 11(5), 387-394. https://doi.org/10.1016/0167-7152(91)90186-U

#<a id="5">[5]</a>
#  Tibshirani, R., Walther, G., & Hastie, T. (2001). Estimating the number of clusters in a data set via the gap statistic. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 63(2), 411-423. https://doi.org/10.1111/1467-9868.00293

#<a id="6">[6]</a>
#  Witten, D. M., & Tibshirani, R. (2010). sparcl: Perform Sparse Hierarchical Clustering and Sparse K-Means Clustering. R package. https://CRAN.R-project.org/package=sparcl

#<a id="7">[7]</a>
#  José Antonio Bernabé-Díaz, Manuel Franco, Juana-María Vivo, Manuel Quesada-Martínez, & Jesualdo T. Fernández-Breis (2022). An automated process for supporting decisions in clustering-based data analysis. Computer Methods and Programs in Biomedicine, 219, 106765. https://doi.org/10.1016/j.cmpb.2022.106765



