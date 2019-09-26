#' @title Stability index.
#' @name stability
#' @aliases stability
#' @description
#' This analysis permits to estimate whether the clustering is meaningfully
#' affected by small variations in the sample. First, a clustering using the
#' k-means algorithm is carried out. The value of \code{k} can be provided by the user.
#' Then, the stability index is the mean of the Jaccard coefficient
#' values of a number of \code{bs} bootstrap replicates. The values are in the range [0,1],
#' having the following meaning:
#' \itemize{
#' \item Unstable: [0, 0.60[.
#' \item Doubtful: [0.60, 0.75].
#' \item Stable: ]0.75, 0.85].
#' \item Highly Stable: ]0.85, 1].
#' }
#'
#' @param data A \code{\link{SummarizedExperiment}}.
#' The SummarizedExperiment must contain an assay with the following structure:
#' A valid header with names. The first  column of the header is the ID or name
#' of the instance of the dataset (e.g., ontology, pathway, etc.) on which the
#' metrics are measured.
#' The other columns of the header contains the names of the metrics.
#' The rows contains the measurements of the metrics for each instance in the dataset.
#' @param k Positive integer. Number of clusters between [2,15] range.
#' @param bs Positive integer. Bootstrap value to perform the resampling.
#' @param getImages Boolean. If true, a plot is displayed.
#'
#' @return A \code{\link{SummarizedExperiment}},
#' containing an assay with the stability measurements and means for 1 to k clusters.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' result <- stability(ontMetrics, k=6, getImages=TRUE)
#'
#' @references
#' \insertRef{milligan1996measuring}{evaluomeR}
#'
#' \insertRef{jaccard1901distribution}{evaluomeR}
#'
#'
stability <- function(data, k=5, bs=100, getImages=TRUE) {

  data <- as.data.frame(assay(data))

  checkKValue(k)

  suppressWarnings(
    runStabilityIndex(data, k.min=k, k.max=k, bs))
  stabilityDataFrame <- suppressWarnings(
    runStabilityIndexTableRange(k.min=k, k.max=k))
  if (getImages == TRUE) {
    suppressWarnings(
      runStabilityIndexK_IMG(bs, k.min = k, k.max = k))
  }
  se <- createSE(stabilityDataFrame)
  return(se)
}


#' @title Stability index for a range of k clusters.
#' @name stabilityRange
#' @aliases stabilityRange
#' @description
#' This analysis permits to estimate whether the clustering is meaningfully
#' affected by small variations in the sample. For a range of k values (\code{k.range}),
#' a clustering using the k-means algorithm is carried out.
#' Then, the stability index is the mean of the Jaccard coefficient
#' values of a number of \code{bs} bootstrap replicates. The values are in the range [0,1],
#' having the following meaning:
#' \itemize{
#' \item Unstable: [0, 0.60[.
#' \item Doubtful: [0.60, 0.75].
#' \item Stable: ]0.75, 0.85].
#' \item Highly Stable: ]0.85, 1].
#' }
#'
#' @inheritParams stability
#' @param k.range Concatenation of two positive integers.
#' The first value \code{k.range[1]} is considered as the lower bound of the range,
#' whilst the second one, \code{k.range[2]}, as the higher. Both values must be
#' contained in [2,15] range.
#'
#' @return A \code{\link{SummarizedExperiment}} containing the stability measurements and
#' means for 1 to k clusters.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' result <- stabilityRange(ontMetrics, k.range=c(2,3))
#'
#' @references
#' \insertRef{milligan1996measuring}{evaluomeR}
#'
#' \insertRef{jaccard1901distribution}{evaluomeR}
#'
#'
stabilityRange <- function(data, k.range=c(2,15), bs=100,
                           getImages=TRUE) {
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
  }

  data <- as.data.frame(assay(data))

  suppressWarnings(
    runStabilityIndex(data, k.min=k.min, k.max=k.max, bs))
  stabilityDataFrame <- suppressWarnings(
    runStabilityIndexTableRange(k.min=k.min, k.max=k.max))

  if (getImages == TRUE) {
    suppressWarnings(
      runStabilityIndexK_IMG(bs, k.min=k.min, k.max=k.max))
    suppressWarnings(
      runStabilityIndexMetric_IMG(bs, k.min=k.min, k.max=k.max))
  }
  se <- createSE(stabilityDataFrame)
  return(se)
}

runStabilityIndex <- function(data, k.min, k.max, bs) {
  inversa=NULL
  m.stab.global = NULL
  todo.estable = NULL
  datos.bruto=data
  names.metr=names(datos.bruto)[-c(1)]

  pkg.env$names.metr = names.metr

  bs.values=c(bs)
  contador=0
  i.min=k.min
  i.max=k.max

  for (i.metr in 1:length(names.metr)) {
    cat("Processing metric: ", names.metr[i.metr],"(", i.metr,")\n")
    m.stab.global[[i.metr]]=matrix(data=NA, nrow=length(names.metr),
                                   ncol=length(i.min:i.max))
    for (j.k in i.min:i.max) {
      cat("\tCalculation of k = ", j.k,"\n")
      estable=NULL
      contador=contador+1
      i=i.metr+1
      estable$n.metric=i.metr
      estable$name.metric=names.metr[i.metr]
      estable$n.k=j.k
      estable$name.ontology=names(datos.bruto[1])

      km5=NULL
      v.size=length(levels(as.factor(datos.bruto[,i])))
      if (v.size>=j.k) {

        km5$cluster=boot.cluster(data=datos.bruto[,i], nk=j.k, B=bs)
        km5$bspart=km5$cluster$partition
        km5$jac=km5$cluster$means

        km5$centr=by(datos.bruto[,i],km5$bspart,mean)
        for (km5.i in 1:length(km5$centr)) {
          km5$means[which(km5$bspart==km5.i)]=km5$centr[km5.i]}

        km5$bspart.or=ordered(km5$means,labels=seq(1,length(km5$centr)))

        km5$bspart.inv=ordered(km5$means,labels=seq(length(km5$centr),1))

        km5$jac.or=km5$jac[order(km5$centr)]
        km5$jac.inv=km5$jac[order(km5$centr,decreasing=TRUE)]

        if (any(inversa==names.metr[i.metr])) {
          km5$partition=km5$bspart.inv
          km5$jac.stab=km5$jac.inv
        } else {
          km5$partition=km5$bspart.or
          km5$jac.stab=km5$jac.or
        }

        m.stab.global[[i.metr]][j.k] = mean(km5$jac.stab)
        estable[[which(bs.values==bs)]]=km5
      } else {
        km5$bspart=rep(NA,length(datos.bruto[,i]))
        km5$jac=rep(NA,j.k)
        km5$centr=rep(NA,j.k)
        km5$means=km5$bspart
        km5$bspart.or=km5$bspart
        km5$bspart.inv=km5$means
        km5$jac.or=km5$jac
        km5$jac.inv=km5$jac
        km5$partition=km5$bspart.inv
        km5$jac.stab=km5$jac.inv
        m.stab.global[[i.metr]][j.k]=mean(km5$jac.stab)
        estable[[which(bs.values==bs)]]=km5
      }
      estable$km5.dynamic=km5$partition
      todo.estable[[contador]]=estable
    }
  }

  e.stab.global=NULL
  for (j.k in i.min:i.max) {
    e.stab.global[[j.k]]=matrix(data=NA, nrow=length(names.metr), ncol=length(i.min:i.max))
    for (i.metr in 1:length(names.metr)) {
      e.stab.global[[j.k]][i.metr]=m.stab.global[[i.metr]][j.k]
    }
  }

  pkg.env$m.stab.global = m.stab.global
  pkg.env$e.stab.global = e.stab.global
  # return(stabilityDataFrame)
  return(NULL)
}

runStabilityIndexTableRange <- function(k.min, k.max) {
  stabilityDataList = list()

  m.stab.global = pkg.env$m.stab.global
  names.metr = pkg.env$names.metr

  # Build header
  header <- list("Metric")
  for (k in k.min:k.max) {
    nextHeader <- paste("Mean_stability_k_", k, sep="")
    header <- c(header, nextHeader)
  }
  header = unlist(header, use.names=FALSE)

  # Build rows
  for (i.metr in 1:length(names.metr)) {
    jackMean = m.stab.global[[i.metr]]
    jackMean = jackMean[!is.na(jackMean)]
    wrapper=NULL # Wrapper object to extract the data
    wrapper$metric=names.metr[i.metr]
    wrapper$jacMean=jackMean
    stabilityDataList[[i.metr]]=unlist(wrapper, use.names=FALSE)
  }  # end for i

  # Transform into dataframe
  stabilityDataFrame = t(data.frame(stabilityDataList))
  colnames(stabilityDataFrame) = header
  rownames(stabilityDataFrame) <- NULL
  return(stabilityDataFrame)
}

# Stability index per K value (x values = matrics)
runStabilityIndexK_IMG <- function(bs, k.min, k.max) {
  ancho=6
  alto=4
  escala=0.6
  escalax=escala
  ajuste=0.5
  escalat=0.5
  escalap=0.4
  contadorFiguras=1
  #Pattern: GlobalStability_K_2, ..., GlobalStability_K_N
  figurename="GlobalStability_K_"

  colores <- c("black", "blue", "red", "darkgreen", "orange", "magenta", "gray")
  ltype <- c(1, 2, 3, 4, 3, 4, 6)

  stype <- c("o", "l")
  pchtype <- c(1, 2, 3, 4, 5, 5)

  e.stab.global = pkg.env$e.stab.global
  names.metr = pkg.env$names.metr
  i.min=k.min
  i.max=k.max
  margins <- par(mar=c(5,5,3,3))
  on.exit(par(margins))
  for (j.k in i.min:i.max) {
    cur.data = e.stab.global[[j.k]]
    cur.data = cur.data[!is.na(cur.data)]
    ymin = min(cur.data)
    xnames=as.character(names.metr)
    ynames="Global Stability Indices"
    g.main=paste(" St. Indices of the metrics for k=", j.k,sep="")
    plot(cur.data, main=g.main, axes=TRUE, col.axis="white",
         xlim=c(0.75,length(xnames)+0.25), xlab="", ylim=c(ymin,1),
         ylab=ynames, col=colores[1],type="o", lwd=1, lty=ltype[1])
    axis(1,at=1:length(xnames),labels=xnames,las=2,cex.axis=0.75)
    axis(2,las=3,cex.axis=0.85)
    labels <- paste("b=", bs, sep = "")
    legend("bottomright", legend=labels, inset=.01, lwd=1, lty=ltype[1:4], col=colores[1:4], cex=0.7, pch=pchtype[1:4])
    mtext(side=1, text="Metrics",line=4)
  }
}

# Stability index per metric (x values = k range)
runStabilityIndexMetric_IMG <- function(bs, k.min, k.max) {
  ancho=6
  alto=4
  escala=0.9
  escalax=escala
  ajuste=0.5
  escalat=0.5
  escalap=0.4
  contadorFiguras=1
  #Pattern: GlobalStability_MetricX, ..., GlobalStability_MetricN
  figurename="GlobalStability_"

  colores <- c("black", "blue", "red", "darkgreen", "orange", "magenta", "gray")
  ltype <- c(1, 2, 3, 4, 3, 4, 6)

  stype <- c("o", "l")
  pchtype <- c(1, 2, 3, 4, 5, 5)

  m.stab.global = pkg.env$m.stab.global
  names.metr = pkg.env$names.metr
  margins <- par(mar=c(5,5,3,3))
  on.exit(par(margins))
  for (i.metr in 1:length(names.metr)) {
    cur.data = m.stab.global[[i.metr]]
    cur.data = cur.data[!is.na(cur.data)]
    ymin = min(cur.data)
    xnames=c(k.min:k.max)
    ynames="Global Stability Indices"
    g.main=paste(" St. Indices of '", names.metr[i.metr], "' for k in [",
                 k.min, ",", k.max,"]",sep="")
    plot(cur.data, main=g.main, axes=TRUE, col.axis="white",
         xlim=c(0.75,length(k.min:k.max)+0.25), xlab="", ylim=c(ymin,1),
         ylab=ynames, col=colores[1],type="o", lwd=1, lty=ltype[1])
    axis(1,at=1:length(k.min:k.max),labels=xnames,las=1,cex.axis=escalax)
    axis(2,las=3,cex.axis=0.85)
    labels <- paste("b=", bs, sep = "")
    legend("bottomright", legend=labels, inset=.01, lwd=1, lty=ltype[1:4], col=colores[1:4], cex=0.7, pch=pchtype[1:4])
    mtext(side=1, text="K values",line=3)
  }

}
