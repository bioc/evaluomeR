% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qualityIndices.R
\name{qualityRange}
\alias{qualityRange}
\title{Goodness of classifications for a range of k clusters.}
\usage{
qualityRange(data, k.range = c(3, 5), cbi = "kmeans",
  getImages = FALSE, all_metrics = FALSE, seed = NULL, ...)
}
\arguments{
\item{data}{A \code{\link{SummarizedExperiment}}.
The SummarizedExperiment must contain an assay with the following structure:
A valid header with names. The first  column of the header is the ID or name
of the instance of the dataset (e.g., ontology, pathway, etc.) on which the
metrics are measured.
The other columns of the header contains the names of the metrics.
The rows contains the measurements of the metrics for each instance in the dataset.}

\item{k.range}{Concatenation of two positive integers.
The first value \code{k.range[1]} is considered as the lower bound of the range,
whilst the second one, \code{k.range[2]}, as the higher. Both values must be
contained in [2,15] range.}

\item{cbi}{Clusterboot interface name (default: "kmeans"):
"kmeans", "clara", "clara_pam", "hclust", "pamk", "pamk_pam", "pamk".
Any CBI appended with '_pam' makes use of \code{\link{pam}}.
The method used in 'hclust' CBI is "ward.D2".}

\item{getImages}{Boolean. If true, a plot is displayed.}

\item{all_metrics}{Boolean. If true, clustering is performed upon all the dataset.}

\item{seed}{Positive integer. A seed for internal bootstrap.}
}
\value{
A list of \code{\link{SummarizedExperiment}} containing the silhouette width measurements and
cluster sizes from \code{k.range[1]} to \code{k.range[2]}. The position on the list matches
with the k-value used in that dataframe. For instance, position 5
represents the dataframe with k = 5.
}
\description{
The goodness of the classifications are assessed by validating the clusters
generated for a range of k values. For this purpose, we use the Silhouette width as validity index.
This index computes and compares the quality of the clustering outputs found
by the different metrics, thus enabling to measure the goodness of the
classification for both instances and metrics. More precisely, this measurement
provides an assessment of how similar an instance is to other instances from
the same cluster and dissimilar to the rest of clusters. The average on all
the instances quantifies how the instances appropriately are clustered. Kaufman
and Rousseeuw suggested the interpretation of the global Silhouette width score
as the effectiveness of the clustering structure. The values are in the
range [0,1], having the following meaning:

\itemize{
\item There is no substantial clustering structure: [-1, 0.25].
\item The clustering structure is weak and could be artificial: ]0.25, 0.50].
\item There is a reasonable clustering structure: ]0.50, 0.70].
\item A strong clustering structure has been found: ]0.70, 1].
}
}
\examples{
# Using example data from our package
data("ontMetrics")
# Without plotting
dataFrameList = qualityRange(ontMetrics, k.range=c(2,3), getImages = FALSE)

}
\references{
\insertRef{kaufman2009finding}{evaluomeR}
}
