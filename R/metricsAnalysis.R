#' @title Minimum and maximum metric values plot.
#' @name plotMetricsMinMax
#' @aliases plotMetricsMinMax
#' @description
#' It plots the minimum, maximum and standard deviation
#' values of the metrics in a \code{\link{SummarizedExperiment}} object.
#'
#' @inheritParams stability
#'
#' @return Nothing.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' plotMetricsMinMax(ontMetrics)
#'
plotMetricsMinMax <- function(data) {
  data <- as.data.frame(assay(data))
  # Prepare data for plotting
  # Data matrix without descritive column
  matrix = data.matrix(data[,-1])
  maxs = matrixStats::colMaxs(matrix)
  mins = matrixStats::colMins(matrix)
  means = colMeans(matrix)
  sd = matrixStats::colSds(matrix)

  dataStats = matrix(NA, nrow=5, ncol = length(data[,-1]), byrow=TRUE,
                     dimnames = list(c("Metric", "Min","Max","Mean","Sd"),
                                     c(colnames(data[,-1]))))
  dataStats["Metric",] = colnames(data[,-1])
  dataStats["Min",] = mins
  dataStats["Max",] = maxs
  dataStats["Mean",] = means
  dataStats["Sd",] = sd
  dataStats.df = as.data.frame(dataStats)
  dataStats.df.t = as.data.frame(t(dataStats.df))
  # Factor to numeric conversion
  dataStats.df.t$Min = as.numeric(as.character(dataStats.df.t$Min))
  dataStats.df.t$Max = as.numeric(as.character(dataStats.df.t$Max))
  dataStats.df.t$Mean = as.numeric(as.character(dataStats.df.t$Mean))
  dataStats.df.t$Sd = as.numeric(as.character(dataStats.df.t$Sd))

  ## Plotting
  p <- ggplot(dataStats.df.t, aes(x=dataStats.df.t$Metric)) +
    geom_linerange(aes(ymin=dataStats.df.t$Min,ymax=dataStats.df.t$Max),
                   linetype=2,color="#4E84C4") +
    geom_point(aes(y=dataStats.df.t$Min),size=3,color="#00AFBB") +
    geom_point(aes(y=dataStats.df.t$Max),size=3,color="#FC4E07") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          #axis.text.y = element_blank(),
          text = element_text(size=15),
          axis.line = element_line(colour = "black",
                                   size = 1, linetype = "solid")
    ) +
    geom_errorbar(aes(ymin=(dataStats.df.t$Max-dataStats.df.t$Sd),
                      ymax=(dataStats.df.t$Max+dataStats.df.t$Sd)), width=.2,
                  position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=(dataStats.df.t$Min-dataStats.df.t$Sd),
                      ymax=(dataStats.df.t$Min+dataStats.df.t$Sd)), width=.2,
                  position=position_dodge(.9)) +
    scale_y_continuous(breaks=seq(round(min(dataStats.df.t$Min-dataStats.df.t$Sd)), # 10 ticks across min - max range
                                  round(max(dataStats.df.t$Max+dataStats.df.t$Sd)),
                                  round((max(dataStats.df.t$Max)-min(dataStats.df.t$Min)))/10),
                       labels=function(x) sprintf("%.2f", x)) + # Two decimals
    labs(x = "Metrics", y = "Metric value", title = "Min/max/sd values across metrics") +
    guides(fill=TRUE)

  print(p)
}


#' @title Metric values as a boxplot.
#' @name plotMetricsBoxplot
#' @aliases plotMetricsBoxplot
#' @description
#' It plots the value of the metrics in a \code{\link{SummarizedExperiment}}
#' object as a boxplot.
#'
#' @inheritParams stability
#'
#' @return Nothing.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' plotMetricsBoxplot(ontMetrics)
#'
plotMetricsBoxplot <- function(data) {
  data <- as.data.frame(assay(data))
  num_metrics_plot=20
  data.metrics = data[,-1] # Removing Description column

  metrics_length = length(colnames(data.metrics))
  num_iterations = round(metrics_length/num_metrics_plot)
  if (num_iterations > 0) {
    num_iterations = num_iterations - 1
  }
  for (iteration in 0:num_iterations) {
    i = 1
    rangeStart = (iteration*num_metrics_plot)+1
    rangeEnd = rangeStart+num_metrics_plot-1
    if (rangeEnd > metrics_length) {
      rangeEnd = metrics_length
    }
    suppressMessages({
      data.melt = melt(data.metrics[,rangeStart:rangeEnd])
    })
    # Melting 1 variable (e.g: data.metrics[,11:11])
    # won't create $variable column in data.melt.
    if (rangeStart == rangeEnd) {
      metricName = data[rangeStart, "Description"]
      data.melt$variable = rep(metricName, length(data.melt$value))
    }
    p <- ggplot(data.melt, aes(x=data.melt$variable, y=data.melt$value)) +
      geom_boxplot(
        #aes(fill=data.melt$variable), # Colors
        outlier.colour = "black",
        outlier.alpha = 0.7,
        outlier.shape = 21,
        show.legend = FALSE
      ) +
      #scale_y_continuous(limits = quantile(data.melt$value, c(0.1, 0.9))) +
      scale_color_grey() +
      theme_bw() +
      theme(
        text = element_text(size=20),
        axis.text.x = element_text(angle = 90)
      ) +
      labs(x = "Metrics", y="Metric value", fill="Metrics")
    # compute lower and upper whiskers
    #ylim1 = boxplot.stats(data.melt$value)$stats[c(1, 5)]
    # scale y limits based on ylim1
    #p1 = p + coord_cartesian(ylim = ylim1*1.05)
    print(p)
  }
}

#' @title Metric values clustering.
#' @name plotMetricsCluster
#' @aliases plotMetricsCluster
#' @description
#' It clusters the value of the metrics in a \code{\link{SummarizedExperiment}}
#' object as a boxplot.
#'
#' @inheritParams stability
#'
#' @return Nothing.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' plotMetricsCluster(ontMetrics)
#'
plotMetricsCluster <- function(data) {
  data <- as.data.frame(assay(data))
  data.metrics = data[,-1] # Removing Description column
  d <- dist(t(data.metrics), method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward.D2")
  theme_set(theme_bw())
  p <- ggdendrogram(fit, rotate = FALSE, size = 2) + # display dendogram
    theme(
      text = element_text(size=15)
    ) +
    labs(title="Metrics dendrogram")
  print(p)
}

#' @title Metric values as violin plot.
#' @name plotMetricsViolin
#' @aliases plotMetricsViolin
#' @description
#' It plots the value of the metrics in a \code{\link{SummarizedExperiment}}
#' object as a violin plot.
#'
#' @inheritParams stability
#'
#' @return Nothing.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' plotMetricsViolin(ontMetrics)
#'
plotMetricsViolin <- function(data) {
  data <- as.data.frame(assay(data))
  data.metrics = data[,-1] # Removing Description column
  num_metrics_plot=20

  metrics_length = length(colnames(data.metrics))
  num_iterations = round(metrics_length/num_metrics_plot)
  if (num_iterations > 0) {
    num_iterations = num_iterations - 1
  }
  for (iteration in 0:num_iterations) {
      i = 1
      rangeStart = (iteration*num_metrics_plot)+1
      rangeEnd = rangeStart+num_metrics_plot-1
      if (rangeEnd > metrics_length) {
        rangeEnd = metrics_length
      }
    suppressMessages({
      data.melt = melt(data.metrics[,rangeStart:rangeEnd])
    })
    # Melting 1 variable (11:11), won't create $variable column in data.melt.
    if (rangeStart == rangeEnd) {
      metricName = data[rangeStart, "Description"]
      data.melt$variable = rep(metricName, length(data.melt$value))
    }
    p <- ggplot(data.melt, aes(x=data.melt$variable, y=data.melt$value)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1) +
      #scale_y_continuous(limits = quantile(data.melt$value, c(0.1, 0.9))) +
      scale_color_grey() +
      theme_bw() +
      theme(
        text = element_text(size=20),
        axis.text.x = element_text(angle = 90)
      ) +
      labs(x = "Metrics", y="Metric value", fill="Metrics")
    # compute lower and upper whiskers
    #ylim1 = boxplot.stats(data.melt$value)$stats[c(1, 5)]
    # scale y limits based on ylim1
    #p1 = p + coord_cartesian(ylim = ylim1*1.05)
    print(p)
  }
}


#
# It returns true if value is in range (0.5, 0.7]
#
isReasonable <- function(value) {
  return(value > 0.5 && value <= 0.7)
}

getLargestSilWidth <- function(qualityDf, metric, k1, k2) {
  k1KSil = qualityDf[metric, k1]
  k2KSil = qualityDf[metric, k2]
  k = NULL
  if (k1KSil >= k2KSil) {
    k = k1
  } else {
    k = k2
  }
  return(getFormattedK(k))
}

#
# It transform a string 'k_X' into 'X'.
# For instace, input is 'k_4', output is '4'
#
getFormattedK <- function(k) {
  return(gsub("^.*_","", k))
}

#' @title Calculating the optimal value of k.
#' getOptimalKValue
#' @aliases getOptimalKValue
#' @description
#' This method finds the optimal value of K per each metric.
#'
#' @param stabData An output \code{\link{SummarizedExperiment}} from
#' a \code{\link{stabilityRange}} execution.
#'
#' @param qualData An output \code{\link{SummarizedExperiment}} from
#' a \code{\link{qualityRange}} execution.
#'
#' @param k.range A range of K values to limit the scope of the
#' analysis.
#'
#' @return It returns a dataframe following the schema:
#' \code{metric}, \code{optimal_k}.
#'
#' @examples
#' # Using example data from our package
#' data("rnaMetrics")
#' stabilityData <- stabilityRange(data=rnaMetrics, k.range=c(2,4), bs=20, getImages = FALSE)
#' qualityData <- qualityRange(data=rnaMetrics, k.range=c(2,4), getImages = FALSE)
#' kOptTable = getOptimalKValue(stabilityData, qualityData)
#'
getOptimalKValue <- function(stabData, qualData, k.range=NULL) {
  checkStabilityQualityData(stabData, qualData)

  if (!is.null(k.range)) {
    k.range.length = length(k.range)
    if (k.range.length != 2) {
      stop("k.range length must be 2")
    }
    k.min = k.range[1]
    k.max = k.range[2]
    checkKValue(k.min)
    checkKValue(k.max)
    if (k.max < k.min) {
      stop("The first value of k.range cannot be greater than its second value")
    } else if (k.min == k.max) {
      stop("Range start point and end point are equals")
    }
  }

  stabDf = standardizeStabilityData(stabData, k.range)
  qualDf = standardizeQualityData(qualData, k.range)

  metrics = rownames(stabDf)
  STABLE_CLASS = 0.75

  outputTable = as.data.frame(metrics)
  rownames(outputTable) = metrics
  outputTable = outputTable[, -1]
  optimalKs = list()
  stabMaxKs = list() # List of maximum K for the stability of metric X
  stabMaxKsStability = list() # Stability of the current K in stabMaxKs
  stabMaxKsQuality = list() # Quality of the current K in stabMaxKs
  qualMaxKs = list() # List of maximum K for the quality of metric X
  qualMaxKsStability = list() # Stability of the current K in qualMaxKs
  qualMaxKsQuality = list() # Quality of the current K in qualMaxKs

  for (metric in metrics) {
    cat("Processing metric: ", metric, "\n")
    stabMaxK = colnames(stabDf[metric, ])[apply(stabDf[metric, ],1,which.max)] # k1
    stabMaxKFormatted = getFormattedK(stabMaxK)
    stabMaxVal = stabDf[metric, stabMaxK]
    qualMaxK = colnames(qualDf[metric, ])[apply(qualDf[metric, ],1,which.max)] # k2
    qualMaxVal = qualDf[metric, qualMaxK]
    qualMaxKFormatted = getFormattedK(qualMaxK)
    ## Info for output table
    stabMaxKs = append(stabMaxKs, stabMaxKFormatted)
    stabMaxKsStability = append(stabMaxKsStability, stabDf[metric, stabMaxK]);
    stabMaxKsQuality = append(stabMaxKsQuality, qualDf[metric, stabMaxK]);

    qualMaxKs = append(qualMaxKs, qualMaxKFormatted)
    qualMaxKsStability = append(qualMaxKsStability, stabDf[metric, qualMaxK]);
    qualMaxKsQuality = append(qualMaxKsQuality, qualDf[metric, qualMaxK]);
    ##
    # CASE 1
    if (length(stabMaxK) >= 1 && length(qualMaxK) >= 1 && identical(stabMaxK, qualMaxK)) {
      k = stabMaxKFormatted
      cat("\tMaximum stability and quality values matches the same K value: '", k ,"'\n")
      optimalKs = append(optimalKs, k)
      # CASE 2
    } else if ((stabMaxVal >= STABLE_CLASS && stabDf[metric, qualMaxK] >= STABLE_CLASS) ||
               (stabMaxVal < STABLE_CLASS && stabDf[metric, qualMaxK] < STABLE_CLASS)) {
      if (stabMaxVal >= STABLE_CLASS && stabDf[metric, qualMaxK] >= STABLE_CLASS) {
        cat("\tBoth Ks have a stable classification: '",
            stabMaxKFormatted, "', '", qualMaxKFormatted ,"'\n")
      } else {
        cat("\tBoth Ks does not have a stable classification: '",
            stabMaxKFormatted, "', '", qualMaxKFormatted ,"'\n")
      }
      k = getLargestSilWidth(qualDf, metric, stabMaxK, qualMaxK)
      cat("\tUsing '", k, "' since it provides higher silhouette width\n")
      optimalKs = append(optimalKs, k)
      # CASE 3
    } else if (stabMaxVal >= STABLE_CLASS && stabDf[metric, qualMaxK] < STABLE_CLASS) {
      cat("\tStability k '", stabMaxKFormatted, "' is stable but quality k '",
          qualMaxKFormatted,"' is not\n")
      kSil = qualDf[metric, stabMaxK]
      if (isReasonable(kSil)) {
        cat("\tUsing ", stabMaxKFormatted, " since its silhouette width is at least reasonable\n")
        optimalKs = append(optimalKs, stabMaxKFormatted)
      } else {
        k = getLargestSilWidth(qualDf, metric, stabMaxK, qualMaxK)
        cat("\tUsing '", k, "' since it provides higher silhouette width\n")
        optimalKs = append(optimalKs, k)
      }
      # CASE 4
    } else if (stabMaxVal < STABLE_CLASS && stabDf[metric, qualMaxK] >= STABLE_CLASS) {
      cat("\t Quality k '", qualMaxKFormatted, "' is stable but stability k '",
          stabMaxKFormatted,"' is not\n")
      kSil = qualDf[metric, qualMaxK]
      if (isReasonable(kSil)) {
        cat("\tUsing ", qualMaxKFormatted, " since its silhouette width is at least reasonable\n")
        optimalKs = append(optimalKs, stabMaxKFormatted)
      } else {
        k = getLargestSilWidth(qualDf, metric, stabMaxK, qualMaxK)
        cat("\tUsing '", k, "' since it provides higher silhouette width\n")
        optimalKs = append(optimalKs, k)
      }
    } else { # This should not happen but it might come in handy to check errors
      optimalKs = append(optimalKs, -1)
    }
  }

  outputTable["Stability_max_k"] = unlist(stabMaxKs)
  outputTable["Stability_max_k_stab"] = unlist(stabMaxKsStability)
  outputTable["Stability_max_k_qual"] = unlist(stabMaxKsQuality)

  outputTable["Quality_max_k"] = unlist(qualMaxKs)
  outputTable["Quality_max_k_stab"] = unlist(qualMaxKsStability)
  outputTable["Quality_max_k_qual"] = unlist(qualMaxKsQuality)

  outputTable["Global_optimal_k"] = unlist(optimalKs)

  return(outputTable)
}

checkStabilityQualityData <- function(stabData, qualData) {
  stabDf = assay(stabData)
  lengthStabDf = length(colnames(stabDf[,-1]))
  stabRangeStart = gsub("^.*_.*_.*_","", colnames(stabDf[,-1])[1]) # Mean_stability_k_2 -> 2
  stabRangeEnd = gsub("^.*_.*_.*_","", colnames(stabDf[,-1])[lengthStabDf])
  lengthQual = length(qualData)
  namesQual = names(qualData)
  qualRangeStart = getFormattedK(namesQual[1]) # k_2 -> 2
  qualRangeEnd = getFormattedK(namesQual[lengthQual])
  if (stabRangeStart != qualRangeStart || stabRangeEnd != qualRangeEnd) {
    stop("Stability data and quality data have different k ranges")
  }
  stabMetricsList = as.character(stabDf[,"Metric"])
  qualMetricsList = as.character(
    assay(getDataQualityRange(qualData, as.numeric(qualRangeStart)))[,"Metric"]
  )
  if (!identical(stabMetricsList, qualMetricsList)) {
    stop("Stability data and quality data have different metrics")
  }
}

#
# It transforms the output of qualityRange method
# into a dataframe like this:
# (rownames)       k_2       k_3       k_4
# DegFact          0.6171262 0.6278294 0.4882649
# ...
# So that the input of getOptimalKValue has always a
# standardized dataframe to process.
#

standardizeQualityData <- function(qualData, k.range=NULL) {
  lengthQuality = length(qualData)
  qualRangeStart = getFormattedK(names(qualData)[1])
  qualRangeEnd = getFormattedK(names(qualData)[lengthQuality])
  Metric = NULL
  kValues = list()
  for (i in seq(qualRangeStart, qualRangeEnd, 1)) {
    curQual = as.data.frame(assay(getDataQualityRange(qualData, i)))
    if (i == qualRangeStart) {
      Metric = as.character(levels(curQual$Metric))
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

#
# It transforms the output of stabilityRange method
# into a dataframe like this:
# (rownames)       k_2       k_3       k_4
# RIN              0.6171262 0.6278294 0.4882649
# ...
# So that the input of getOptimalKValue has always a
# standardized dataframe to process.
#
standardizeStabilityData <- function(stabData, k.range=NULL) {
  stabDf = as.data.frame(assay(stabData))
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
