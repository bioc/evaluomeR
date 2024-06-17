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
#' object a an hclust dendogram from \code{\link{stats}}. By default distance is measured in 'euclidean'
#' and hclust method is 'ward.D20.
#'
#' @inheritParams stability
#' @param scale Boolean. If true input data is scaled. Default: FALSE.
#' @param k Integer. If not NULL a 'cutree' cut on the cluster is done. Default: NULL
#'
#' @return An hclust object.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' plotMetricsCluster(ontMetrics, scale=TRUE)
#'
plotMetricsCluster <- function(data, scale=FALSE, k=NULL) {
  data <- as.data.frame(assay(data))
  data.metrics = data[,-1] # Removing Description column
  if (isTRUE(scale)) {
    data.metrics = base::scale(data.metrics)
  }
  d <- dist(t(data.metrics), method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward.D2")

  dend <- as.dendrogram(fit)

  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), cex = 0.7, col = "blue")
  plot(dend, xlab = "", sub="", ylab = "Euclidean distance",
       main = paste0("Metrics dendrogram 'ward.D2'"), nodePar = nodePar)

  if (!is.null(k)) {
    dendextend::rect.dendrogram(dend , k=k, border="purple")
  }
  return(fit)
}

#' @title Metric values as violin plot.
#' @name plotMetricsViolin
#' @aliases plotMetricsViolin
#' @description
#' It plots the value of the metrics in a \code{\link{SummarizedExperiment}}
#' object as a violin plot.
#'
#' @inheritParams stability
#' @param nplots Positive integer. Number of metrics per violin plot. Default: 20.
#'
#' @return Nothing.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' plotMetricsViolin(ontMetrics)
#'
plotMetricsViolin <- function(data, nplots=20) {
  data <- as.data.frame(assay(data))
  data.metrics = data[,-1] # Removing Description column
  nplots=20

  metrics_length = length(colnames(data.metrics))
  num_iterations = round(metrics_length/nplots)
  if (num_iterations > 0) {
    num_iterations = num_iterations - 1
  }
  for (iteration in 0:num_iterations) {
      i = 1
      rangeStart = (iteration*nplots)+1
      rangeEnd = rangeStart+nplots-1
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
#' @param stabData An output \code{\link{ExperimentList}} from
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

  metrics = as.character(as.data.frame(assay(stabData))$Metric)
  STABLE_CLASS = 0.75

  outputTable = as.data.frame(metrics)
  #rownames(outputTable) = metrics
  outputTable = outputTable[, -1]
  optimalKs = list()
  stabMaxKs = list() # List of maximum K for the stability of metric X
  stabMaxKsStability = list() # Stability of the current K in stabMaxKs
  stabMaxKsQuality = list() # Quality of the current K in stabMaxKs
  qualMaxKs = list() # List of maximum K for the quality of metric X
  qualMaxKsStability = list() # Stability of the current K in qualMaxKs
  qualMaxKsQuality = list() # Quality of the current K in qualMaxKs

  for (metric in metrics) {
    message("Processing metric: ", metric, "\n")
    stabMaxK = colnames(stabDf[metric, ])[apply(stabDf[metric, ],1,which.max)] # ks
    stabMaxKFormatted = getFormattedK(stabMaxK)
    stabMaxVal = stabDf[metric, stabMaxK]
    qualMaxK = colnames(qualDf[metric, ])[apply(qualDf[metric, ],1,which.max)] # kg
    qualMaxKFormatted = getFormattedK(qualMaxK)
    qualMaxVal = qualDf[metric, qualMaxK]
    ## Info for output table
    stabMaxKs = append(stabMaxKs, stabMaxKFormatted)
    stabMaxKsStability = append(stabMaxKsStability, stabDf[metric, stabMaxK]);
    stabMaxKsQuality = append(stabMaxKsQuality, qualDf[metric, stabMaxK]);

    qualMaxKs = append(qualMaxKs, qualMaxKFormatted)
    qualMaxKsStability = append(qualMaxKsStability, stabDf[metric, qualMaxK]);
    qualMaxKsQuality = append(qualMaxKsQuality, qualDf[metric, qualMaxK]);

    # CASE 1: ks == kg
    if (identical(stabMaxK, qualMaxK)) {
      k = stabMaxKFormatted
      message("\tMaximum stability and quality values matches the same K value: '", k ,"'\n")
      optimalKs = append(optimalKs, k)
    } else {
      # CASE 2: ks != kg
      if (stabMaxVal > STABLE_CLASS && stabDf[metric, qualMaxK] > STABLE_CLASS) {
        # Both stables
        message("\tBoth Ks have a stable classification: '",
            stabMaxKFormatted, "', '", qualMaxKFormatted ,"'\n")
        k = qualMaxKFormatted
        optimalKs = append(optimalKs, k)
        message("\tUsing '", k, "' since it provides higher silhouette width\n")
      } else {
        if (stabMaxVal <= STABLE_CLASS && stabDf[metric, qualMaxK] <= STABLE_CLASS) {
          # Both not stables: S_ks <= 0.75 && S_kg <= 0.75
          message("\tBoth Ks do not have a stable classification: '",
              stabMaxKFormatted, "', '", qualMaxKFormatted ,"'\n")
          k = qualMaxKFormatted
          optimalKs = append(optimalKs, k)
          message("\tUsing '", k, "' since it provides higher silhouette width\n")
        } else {
          # S_ks > 0.75 && Sil_ks > 0.5 && S_kg <= 0.75
          if ((stabMaxVal > STABLE_CLASS) && (qualDf[metric, stabMaxK] > 0.5)
              && (stabDf[metric, qualMaxK] <= STABLE_CLASS)) {
            message("\tStability k '", stabMaxKFormatted, "' is stable but quality k '",
                qualMaxKFormatted,"' is not\n")
            k = stabMaxKFormatted
            optimalKs = append(optimalKs, k)
            message("\tUsing '", k, "' since it provides higher stability\n")
          } else {
            # CASE 3
            if (stabMaxVal > STABLE_CLASS && qualDf[metric, stabMaxK] <= 0.5
                && stabDf[metric, qualMaxK] <= STABLE_CLASS)  {
              message("\tStability k '", stabMaxKFormatted, "' is stable but its silhouette value is not reasonable\n")
              if (qualMaxVal > 0.5) { # S_kg > 0.5
                k = qualMaxKFormatted
                optimalKs = append(optimalKs, k)
                message("\tUsing quality '", k, "' since its at least reasonable\n")
              } else {# S_kg <= 0.5
                k = stabMaxKFormatted
                optimalKs = append(optimalKs, k)
                message("\tUsing stability '", k, "' since quality k is not reasonable\n")
              }
            } else { # This should not happen but it might come in handy to check errors
              message("\tUnknown case\n")
              optimalKs = append(optimalKs, -1)
            }
          }
        }
      }
    }
  }

  outputTable["Metric"] = metrics
  outputTable["Stability_max_k"] = unlist(stabMaxKs)
  outputTable["Stability_max_k_stab"] = unlist(stabMaxKsStability)
  outputTable["Stability_max_k_qual"] = unlist(stabMaxKsQuality)

  outputTable["Quality_max_k"] = unlist(qualMaxKs)
  outputTable["Quality_max_k_stab"] = unlist(qualMaxKsStability)
  outputTable["Quality_max_k_qual"] = unlist(qualMaxKsQuality)

  outputTable["Global_optimal_k"] = unlist(optimalKs)

  return(outputTable)

}


#' @title Comparison between two clusterings as plot.
#' plotMetricsClusterComparison
#' @aliases plotMetricsClusterComparison
#' @description
#' It plots a clustering comparison between two different
#' k-cluster vectors for a set of metrics.
#'
#' @inheritParams stability
#' @param k.vector1 Vector of positive integers representing \code{k} clusters.
#' The \code{k} values must be contained in [2,15] range.
#' @param k.vector2 Optional. Vector of positive integers representing \code{k} clusters.
#' The \code{k} values must be contained in [2,15] range.
#'
#' @return Nothing.
#'
#' @examples
#' # Using example data from our package
#' data("rnaMetrics")
#' stabilityData <- stabilityRange(data=rnaMetrics, k.range=c(2,4), bs=20, getImages = FALSE)
#' qualityData <- qualityRange(data=rnaMetrics, k.range=c(2,4), getImages = FALSE)
#' kOptTable = getOptimalKValue(stabilityData, qualityData)
#'
#'
plotMetricsClusterComparison <- function(data, k.vector1, k.vector2=NULL, seed=NULL) {
  if (is.null(seed)) {
    seed = pkg.env$seed
  }
  if (identical(k.vector1, k.vector2)) {
    stop("k.vector1 and k.vector2 are identical")
  }

  data <- as.data.frame(SummarizedExperiment::assay(data))

  numMetrics = length(colnames(data))-1

  if (length(k.vector1) == 1) {
    k.vector1=rep(k.vector1, numMetrics)
  }

  if (is.null(k.vector2)) {
    k.vector2 = k.vector1 # This will colour elipses around the same clusters of k.vector1
  }

  if (length(k.vector2) == 1) {
    k.vector2=rep(k.vector2, numMetrics)
  }

  if (numMetrics != length(k.vector1) || numMetrics != length(k.vector2)
      || length(k.vector1) != length(k.vector2)) {
    stop("Input parameters have different lengths")
  }
  for (i in 1:length(k.vector1)) {
    checkKValue(k.vector1[i])
    checkKValue(k.vector2[i])
  }

  data.metrics=NULL; names.metr=NULL; names.index=NULL;
  k.cl=NULL; k.min=NULL; k.max=NULL;
  data.metrics=NULL; datos.csv=NULL; datos.raw=NULL;
  ranges=NULL; mins=NULL; data.l=NULL; data.ms=NULL; k.sig=NULL; k.op.sig=NULL;

  datos.csv = data
  data.metrics <- datos.csv[,-1]
  names.metr <- colnames(datos.csv[,-1])  #nombres de metricas
  names.ont <- datos.csv[,1]

  ranges <- apply(data.metrics, 2, sample.range)
  mins <- apply(data.metrics, 2, sample.min)
  data.l <- sweep(data.metrics, 2, mins, FUN="-")
  data.ms <- sweep(data.l, 2, ranges, FUN="/")

  kcolors=c("black","red","blue","green","magenta","pink","yellow","orange","brown","cyan","gray","darkgreen")

  par(mar=c(4,6,3,3))
  plot(0,0, xlim=range(data.ms), ylim=c(0,length(names.metr)+1),
       lwd=NULL, xlab="", ylab="", xaxt="n", yaxt="n", type="n")
  axis(side=2, at = seq(1,length(names.metr)), labels=names.metr, las=2, cex.axis=.7)
  title(xlab=paste("Scaled raw scores", sep=""), line=1)
  title(ylab="Metrics", line=5)

  for (i.metr in 1:length(names.metr)) { # i.metr= n de metrica #ejemplo
    #  i.metr=1

    i=NULL; clusterk5=NULL; clusterkopt=NULL;
    k.cl=NULL; k.op=NULL; data.plot=NULL;

    i=i.metr

    #kmeans with k.cl classes
    k.cl=k.vector2[i]
    set.seed(seed)
    clusterk5=kmeans(data.ms[,i], centers=k.cl, iter.max = 100)
    ##
    clusterk5$means=by(data.ms[,i],clusterk5$cluster,mean) #calcula las k medias (centroides)
    for (i.5 in 1:length(clusterk5$means)) {
      clusterk5$partition[which(clusterk5$cluster==i.5)]=clusterk5$centers[i.5]
      #asigna valor centroide a todo miembro del cluster
    }
    #Ordenacion de la particion segun el sentido de la metrica (directa/inversa)
    clusterk5$ordered=ordered(clusterk5$partition,labels=seq(1,length(clusterk5$centers)))
    clusterk5$ordered.inv=ordered(clusterk5$partition,labels=seq(length(clusterk5$centers),1))
    clusterk5$partition=clusterk5$ordered
    clusterk5$means=sort(clusterk5$means,decreasing=FALSE)

    #kmeans with k.op classes
    k.op=k.vector1[i]
    set.seed(seed)
    clusterkopt=kmeans(data.ms[,i], centers=k.op, iter.max = 100)

    clusterkopt$means=by(data.ms[,i],clusterkopt$cluster,mean) #calcula las k medias (centroides)
    for (i.opt in 1:length(clusterkopt$means)) {
      clusterkopt$partition[which(clusterkopt$cluster==i.opt)]=clusterkopt$centers[i.opt]
      #asigna valor centroide a todo el cluster
    }

    clusterkopt$ordered=ordered(clusterkopt$partition,labels=seq(1,length(clusterkopt$centers)))
    clusterkopt$ordered.inv=ordered(clusterkopt$partition,labels=seq(length(clusterkopt$centers),1))
    clusterkopt$partition=clusterkopt$ordered
    clusterkopt$means=sort(clusterkopt$means,decreasing=FALSE)

    data.plot=data.frame(data.ms[,i],clusterk5$partition,clusterkopt$partition)
    colnames(data.plot)=c(names.metr[i],"k=5","k_op")
    rownames(data.plot)=names.ont

    xi=data.plot[[1]]
    yi=rep(i.metr,length(xi))
    ci=data.plot[[2]]
    ci=levels(ci)[ci]
    points(xi,yi,type="p", col=kcolors[as.numeric(ci)],lty=1, lwd=1)

    cj=data.plot[[3]]
    for (ellip.j in unique(cj)) {
      xj=mean(range(xi[which(cj==ellip.j)])) #clusterk5$means[ellip.j]
      yj=rep(i.metr,length(xj))
      aj=diff(range(xi[which(cj==ellip.j)]))/2
      draw.ellipse(x=xj, y=yj, a=aj, b=0.3, nv=100,
                   border=kcolors[as.numeric(ellip.j)], lty=1, lwd=2)
    }

  } #end for i.metr
}

checkStabilityQualityData <- function(stabData, qualData) {
  stabDf = assay(stabData) # Getting first assay, which is 'stabData$stability_mean'
  lengthStabDf = length(colnames(stabDf)[-1])
  stabRangeStart = gsub("^.*_.*_.*_","", colnames(stabDf)[-1][1]) # Mean_stability_k_2 -> 2
  stabRangeEnd = gsub("^.*_.*_.*_","", colnames(stabDf)[-1][lengthStabDf])
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
    if (length(values) < length(Metric)) {
      for (i in seq(length(values), length(Metric)-1,1)) {
        values = append(values, NaN)
      }
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

#' @title Calculate the cluster ID from the optimal cluster per metric for each individual.
#' annotateClustersByMetric
#' @aliases annotateClustersByMetric
#' @description
#' Return a named list, where each metric name is linked to a data frame
#' containing the evaluated individuals, their score for the specified metric,
#' and the cluster id in which each individual is classified. This cluster
#' assignment is performed by calculating the optimal k value by evaluome.
#'
#' @param df Input data frame. The first column denotes the identifier of the
#' evaluated individuals. The remaining columns contain the metrics used to
#' evaluate the individuals. Rows with NA values will be ignored.
#' @param k.range Range of k values in which the optimal k will be searched
#' @param bs Bootstrap re-sample param.
#' @param seed Random seed to be used.
#'
#' @return A named list resulting from computing the optimal cluster for each
#' metric. Each metric is a name in the named list, and its content is a
#' data frame that includes the individuals, the value for the corresponding
#' metric, and the cluster id in which the individual has been asigned according
#' to the optimal cluster.
#' @export
#'
#' @examples
#' data("ontMetrics")
#' annotated_clusters=annotateClustersByMetric(ontMetrics, k.range=c(2,3), bs=20, seed=100)
#' View(annotated_clusters[['ANOnto']])
annotateClustersByMetric <- function(df, k.range, bs, seed){
  if (is.null(seed)) {
    seed = pkg.env$seed
  }
  df <- as.data.frame(assay(df))
  # Create a dataframe by removing NAs from the original data.
  df_clean = na.omit(df)

  # Compute stability, quality and the optimal k value from evaluome
  stabilityData <- stabilityRange(data=df_clean, k.range=k.range,
                                  bs=bs, getImages = FALSE, seed=seed)

  qualityData <- qualityRange(data=df_clean, k.range=k.range,
                              getImages = FALSE, seed=seed)

  kOptTable <- getOptimalKValue(stabilityData, qualityData)

  # Get the clusters obtained by evaluome for each k, together with
  # the optimal k
  clusters = as.data.frame(assay(stabilityData$cluster_partition))
  clusters$optimal_k = kOptTable$Global_optimal_k

  # Compose the results
  # For each metric, get the optimal k, get the clusters formed by using
  # that optimal k, include this information in a dataframe
  result_list = list()

  for (i in 1:nrow(clusters)){
    metric = as.character(clusters$Metric[i])
    # Get the optimal k
    optimal_k = clusters$optimal_k[i]

    # Get the clusters formed by using the optimal k as a vector of integers
    optimal_cluster = dplyr::select(clusters, dplyr::contains(as.character(optimal_k)))
    optimal_cluster = optimal_cluster[i,1]
    optimal_cluster = as.character(optimal_cluster)
    optimal_cluster = as.numeric(strsplit(optimal_cluster, ", ")[[1]])

    # Create a dataframe including the individual id, the concerning metric
    # and the cluster id in which the individual is classfied.
    annotated_df_clean = dplyr::select(df_clean, 1)
    annotated_df_clean$cluster = optimal_cluster

    # Merge this dataframe with the original one, so that original individuals
    # removed due to NAs will be present an NA as cluster.
    # Include this dataframe in the named list, using the name of the metric as
    # a key.
    result_list[[metric]] = merge(dplyr::select(df, 1, dplyr::contains(metric)), annotated_df_clean, all.x = TRUE)
  }
  return(result_list)
}


#' @title Get the range of each metric per cluster from the optimal cluster.
#' getMetricRangeByCluster
#' @aliases getMetricRangeByCluster
#' @description
#' Obtains the ranges of the metrics obtained by each optimal cluster.
#'
#' @param df Input data frame. The first column denotes the identifier of the
#' evaluated individuals. The remaining columns contain the metrics used to
#' evaluate the individuals. Rows with NA values will be ignored.
#' @param k.range Range of k values in which the optimal k will be searched
#' @param bs Bootstrap re-sample param.
#' @param seed Random seed to be used.
#'
#' @return A dataframe including the min and the max value for each
#' pair (metric, cluster).
#' @export
#'
#' @examples
#' data("ontMetrics")
#' ranges = getMetricRangeByCluster(ontMetrics, k.range=c(2,3), bs=20, seed=100)
#' View(ranges)
getMetricRangeByCluster <- function(df, k.range, bs, seed) {
  if (is.null(seed)) {
    seed = pkg.env$seed
  }
  df <- as.data.frame(assay(df))
  annotated_clusters_by_metric = annotateClustersByMetric(df, k.range, bs, seed)

  metrics = c()
  cluster_ids = c()
  min_values = c()
  max_values = c()
  # For each metric
  for (metric in names(annotated_clusters_by_metric)){
    # Get the optimal clusters
    annotated_clusters = annotated_clusters_by_metric[[metric]]

    # For each cluster, get the minimal and the maximal value
    for (cluster_id in 1:max(annotated_clusters$cluster, na.rm=T)) {
      concrete_cluster_values = dplyr::filter(annotated_clusters, cluster==cluster_id) %>% dplyr::pull(metric)
      metrics = c(metrics, metric)
      cluster_ids = c(cluster_ids, cluster_id)
      min_values = c(min_values, min(concrete_cluster_values))
      max_values = c(max_values, max(concrete_cluster_values))
    }
  }
  return(data.frame(metric=metrics, cluster=cluster_ids, min_value=min_values, max_value=max_values))
}


#' @title Get the range of each metric per cluster from the optimal cluster.
#' getMetricRangeByCluster
#' @aliases getMetricRangeByCluster
#' @description
#' Obtains the ranges of the metrics obtained by each optimal cluster.
#'
#' @param df Input data frame. The first column denotes the identifier of the
#' evaluated individuals. The remaining columns contain the metrics used to
#' evaluate the individuals. Rows with NA values will be ignored.
#' @param k K value (number of clusters)
#' @param alpha 0 <= alpha <= 1, the proportion of the cases to be trimmed in robust sparse K-means, see \code{\link{RSKC}}.
#' @param L1 A single L1 bound on weights (the feature weights), see \code{\link{RSKC}}.
#' @param seed Random seed to be used.
#'
#' @return A dataframe including the min and the max value for each
#' pair (metric, cluster).
#' @export
#'
#' @examples
#' data("ontMetrics")
#' metricsRelevancy = getMetricsRelevancy(ontMetrics, k=3, alpha=0.1, seed=100)
#' metricsRelevancy$rskc # RSKC output object
#' metricsRelevancy$trimmed_cases # Trimmed cases from input (row indexes)
#' metricsRelevancy$relevancy # Metrics relevancy table
#'
getMetricsRelevancy <- function(df, k, alpha=NULL, L1=NULL, seed=NULL) {
    if (is.null(seed)) {
    seed = pkg.env$seed
  }
  if (is.null(alpha)) {
    alpha = 0.1
  }

  df <- as.data.frame(assay(df))

  if (is.null(L1)) {
    print(paste0("No L1 provided. Computing best L1 boundry with 'sparcl::KMeansSparseCluster.permute'"))
    dataMatrix = as.matrix(df)
    wbounds = seq(2,sqrt(ncol(dataMatrix)), len=30)
    km.perm <- sparcl::KMeansSparseCluster.permute(dataMatrix,K=k,wbounds=wbounds,nperms=5,silent=TRUE)
    L1 = km.perm$bestw

  }


  # Compute RSKC
  print(paste0("Alpha set as: ", alpha))
  print(paste0("L1 set as: ", L1))
  rskc_out = RSKC(df, k, alpha, L1 = L1, nstart = 200,
                  silent=TRUE, scaling = FALSE, correlation = FALSE)
  # Get trimmed cases from input
  union_vector = c(rskc_out$oE,rskc_out$oW)
  union_vector_unique = unique(union_vector)
  union_vector_unique = sort(union_vector_unique)

  # Metrics relevancy
  columns = c('metric', 'weight')
  rskc_df = data.frame(matrix(ncol = length(columns), nrow = length(rskc_out$weights)))
  colnames(rskc_df) = columns
  rskc_df['metric'] = names(rskc_out$weights)
  rskc_df['weight'] = rskc_out$weights
  rskc_df_sorted = rskc_df[order(rskc_df$weight, decreasing = TRUE), ] # Sorting from greater values to lower



  output = NULL
  output$rskc = rskc_out
  output$trimmed_cases = union_vector_unique
  output$relevancy = rskc_df_sorted
  return (output)
}
