#' @title Goodness of classifications.
#' @name quality
#' @aliases quality
#' @description
#' The goodness of the classifications are assessed by validating the clusters
#' generated. For this purpose, we use the Silhouette width as validity index.
#' This index computes and compares the quality of the clustering outputs found
#' by the different metrics, thus enabling to measure the goodness of the
#' classification for both instances and metrics. More precisely, this goodness measurement
#' provides an assessment of how similar an instance is to other instances from
#' the same cluster and dissimilar to all the other clusters. The average on all
#' the instances quantifies how appropriately the instances are clustered. Kaufman
#' and Rousseeuw suggested the interpretation of the global Silhouette width score
#' as the effectiveness of the clustering structure. The values are in the
#' range [0,1], having the following meaning:
#'
#' \itemize{
#' \item There is no substantial clustering structure: [-1, 0.25].
#' \item The clustering structure is weak and could be artificial: ]0.25, 0.50].
#' \item There is a reasonable clustering structure: ]0.50, 0.70].
#' \item A strong clustering structure has been found: ]0.70, 1].
#' }
#'
#' @inheritParams stability
#'
#' @return A \code{\link{SummarizedExperiment}} containing the silhouette width measurements and
#' cluster sizes for cluster \code{k}.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' result = quality(ontMetrics, k=4)
#'
#' @references
#' \insertRef{kaufman2009finding}{evaluomeR}
#'
quality <- function(data, k=5, getImages=TRUE) {

  checkKValue(k)

  data <- getAssay(data, 1)

  suppressWarnings(
    runQualityIndicesSilhouette(data, k.min = k, k.max = k, bs = 1))
  silhouetteDataFrame = suppressWarnings(
    runSilhouetteTable(data, k = k))
  if (getImages == TRUE) {
    suppressWarnings(
      runQualityIndicesSilhouetteK_IMG(k.min = k, k.max = k))
    suppressWarnings(
      runSilhouetteIMG(data, k))
  }
  se <- createSE(silhouetteDataFrame)
  return(se)

}

#' @title Goodness of classifications for a range of k clusters.
#' @name qualityRange
#' @aliases qualityRange
#' @description
#' The goodness of the classifications are assessed by validating the clusters
#' generated for a range of k values. For this purpose, we use the Silhouette width as validity index.
#' This index computes and compares the quality of the clustering outputs found
#' by the different metrics, thus enabling to measure the goodness of the
#' classification for both instances and metrics. More precisely, this measurement
#' provides an assessment of how similar an instance is to other instances from
#' the same cluster and dissimilar to the rest of clusters. The average on all
#' the instances quantifies how the instances appropriately are clustered. Kaufman
#' and Rousseeuw suggested the interpretation of the global Silhouette width score
#' as the effectiveness of the clustering structure. The values are in the
#' range [0,1], having the following meaning:
#'
#' \itemize{
#' \item There is no substantial clustering structure: [-1, 0.25].
#' \item The clustering structure is weak and could be artificial: ]0.25, 0.50].
#' \item There is a reasonable clustering structure: ]0.50, 0.70].
#' \item A strong clustering structure has been found: ]0.70, 1].
#' }
#'
#' @inheritParams stability
#' @param k.range Concatenation of two positive integers.
#' The first value \code{k.range[1]} is considered as the lower bound of the range,
#' whilst the second one, \code{k.range[2]}, as the higher. Both values must be
#' contained in [2,15] range.
#'
#' @return A list of \code{\link{SummarizedExperiment}} containing the silhouette width measurements and
#' cluster sizes from \code{k.range[1]} to \code{k.range[2]}. The position on the list matches
#' with the k-value used in that dataframe. For instance, position 5
#' represents the dataframe with k = 5.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' # Without plotting
#' dataFrameList = qualityRange(ontMetrics, k.range=c(2,6), getImages = FALSE)
#'
#' @references
#' \insertRef{kaufman2009finding}{evaluomeR}
#'
qualityRange <- function(data, k.range=c(3,5), getImages=TRUE) {

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

  data <- getAssay(data, 1)

  suppressWarnings(
    runQualityIndicesSilhouette(data, k.min = k.min, k.max = k.max, bs = 1))
  silhouetteData =  suppressWarnings(
    runSilhouetteTableRange(data, k.min = k.min, k.max = k.max))

  if (getImages == TRUE) {
    suppressWarnings(
      runQualityIndicesSilhouetteK_IMG(k.min = k.min, k.max = k.max))
    suppressWarnings(
      runQualityIndicesSilhouetteMetric_IMG(k.min = k.min, k.max = k.max))
  }
  seList <- createSEList(silhouetteData)
  return(seList)
}

runQualityIndicesSilhouette <- function(data, k.min, k.max, bs) {
  datos.bruto=NULL
  names.metr=NULL
  names.index=NULL

  datos.bruto=data
  names.metr=names(datos.bruto)[-c(1)]
  pkg.env$names.metr = names.metr
  names.index=c("sil")
  pkg.env$names.index = names.index
  k.min=k.min
  k.max=k.max


  estable=NULL
  m.global=NULL
  e.global=NULL
  i.min=NULL
  i.max=NULL
  contador=NULL
  contador=0
  remuestreo=bs

  i.min=k.min
  i.max=k.max

  for (i.metr in 1:length(names.metr)) {
    cat("Processing metric: ", names.metr[i.metr],"(", i.metr,")\n")
    m.global[[i.metr]]=matrix(data=NA, nrow=i.max, ncol=length(i.min:i.max))
    for (j.k in i.min:i.max) {
      cat("\tCalculation of k = ", j.k,"\n")
      e.res=NULL
      e.res.or=NULL
      i=NULL
      j=NULL
      contador=contador+1
      i=i.metr+1

      j=j.k
      #e.res=ClusterStability(dat=datos.bruto[,i],k=j, replicate=remuestreo, type='kmeans')

      e.res$n=contador
      e.res$n.metric=i.metr
      e.res$name.metric=names.metr[i.metr]

      e.res$n.k=j.k
      e.res$name.ontology=datos.bruto$Description

      e.res$km5.dynamic.bs <- quiet(clusterboot(datos.bruto[,i], B=remuestreo, bootmethod="boot",clustermethod=kmeansCBI,krange=5)$partition)
      e.res$kmk.dynamic.bs <- quiet(clusterboot(datos.bruto[,i], B=remuestreo, bootmethod="boot",clustermethod=kmeansCBI,krange=j)$partition)

      e.res.or$centr=by(datos.bruto[,i],e.res$kmk.dynamic.bs,mean)

      for (e.res.or.i in 1:length(e.res.or$centr)) {
        e.res.or$means[which(e.res$kmk.dynamic.bs==e.res.or.i)]=e.res.or$centr[e.res.or.i]}

      e.res$kmk.dynamic.bs.or=ordered(e.res.or$means,labels=seq(1,length(e.res.or$centr)))

      ## Using Silhouette width as index
      metric.onto=datos.bruto[,i.metr+1]
      part.onto=as.numeric(e.res$kmk.dynamic.bs.or)
      sil.w=silhouette(part.onto, dist(metric.onto))
      sil.c=cluster.stats(dist(metric.onto), part.onto, noisecluster=TRUE)
      e.res$sil.w = sil.w
      e.res$sil.c = sil.c
      estable[[contador]]=e.res

      m.global[[i.metr]][j.k,] = mean(sil.w[,"sil_width"])
      #m.global[[i.metr]][j.k,] = e.res$ST_global_sil
    }
  }
  for (j.k in i.min:i.max) {
    e.global[[j.k]]=matrix(data=NA, nrow=length(names.metr), ncol=length(i.min:i.max))
    for (i.metr in 1:length(names.metr)) {
      e.global[[j.k]][i.metr,]=m.global[[i.metr]][j.k,]
    }
  }

  pkg.env$m.global = m.global
  pkg.env$e.global = e.global
  pkg.env$estable = estable
}

# Silhouette width per k (x values = metrics)
runQualityIndicesSilhouetteK_IMG <- function(k.min, k.max) {
  ancho=NULL
  alto=NULL
  ajuste=NULL
  escala=NULL
  escalax=NULL
  escalal=NULL
  escalat=NULL
  escalap=NULL
  ancho=6
  alto=4
  escala=0.75
  escalax=escala
  escalal=0.85
  ajuste=0.5
  escalat=0.5
  escalap=0.4
  listaFigurasGrafMtr=NULL
  contadorFiguras=1

  e.mat.global=NULL
  nombres=NULL
  i.min=NULL
  i.max=NULL
  x=NULL
  x.label=NULL
  x.name=NULL
  y.label=NULL
  m.global = pkg.env$m.global
  e.global = pkg.env$e.global
  e.mat.global=e.global

  names.index = pkg.env$names.index
  i.min=1
  i.max=k.max-(k.min-1)

  names.metr = pkg.env$names.metr
  x=seq(1,length(names.metr))
  x.label="Metrics"
  x.name=names.metr
  y.label="Silhouette avg. width"
  #Pattern: QualityIndices_K_2, ..., QualityIndices_K_N
  figurename="QualityIndices_K_"

  i.min=k.min
  i.max=k.max
  margins <- par(mar=c(5,5,3,3))
  on.exit(par(margins))
  for (m.g in i.min:i.max) {
    g.main=NULL
    leg.g=NULL
    ymin=NULL
    ymax=NULL
    xmin=NULL
    xmax=NULL
    xleg=NULL
    yleg=NULL
    c.max=NULL
    t.linea=NULL
    t.color=NULL
    ymarcas=NULL
    k.classes=NULL
    xmin=min(x)-0.25
    xmax=max(x)+0.25
    xleg=((xmax-xmin)*escalal)+3.2
    c.max=dim(e.mat.global[[m.g]])[2]
    ymin=min(e.mat.global[[m.g]])
    ymax=1
    ymarcas=round(seq(ymin,ymax,length.out=5),2)
    yleg=ymin+((ymax-ymin)/2)*seq(c.max,1,-1)/(2*c.max)
    t.linea=seq(1,c.max)
    t.color=rep("black",c.max)

    k.classes=m.g
    g.main=paste(" Qual. Indices of the metrics for k=",k.classes,sep="")

    par(new=FALSE,bg="white",fg="black")

    for (m in length(names.index)) {
      y=NULL
      y.name=NULL
      y=e.mat.global[[m.g]][,m]
      y.name=names.index[m]
      leg.g[m] <- paste(y.name," avg. width",sep="")
      plot(x,y, type="l", xaxt="n", yaxt="n", xlab="", ylab="", main=g.main, xlim=c(xmin,xmax), ylim=c(ymin,ymax), lty=t.linea[m], col=t.color[m])
      par(new=TRUE)
      plot(x,y, type="o", xaxt="n", yaxt="n", xlab="", ylab="", main=g.main, xlim=c(xmin,xmax), ylim=c(ymin,ymax), lty=t.linea[m], col=t.color[m])
      par(new=TRUE)
    }
    mtext(side=1, text=x.label,line=4)
    mtext(side=2, text=y.label,line=3)
    axis(side=1, at=x, labels=x.name, las=3, cex.axis=escalax)
    axis(side=2, at=ymarcas, labels=ymarcas, cex.axis=escalal)
    par(new=FALSE)
  }
}

# Silhouette width per metric (x values = k range)
runQualityIndicesSilhouetteMetric_IMG <- function(k.min, k.max) {
  ancho=NULL
  alto=NULL
  ajuste=NULL
  escala=NULL
  escalax=NULL
  escalal=NULL
  escalat=NULL
  escalap=NULL
  ancho=6
  alto=4
  escala=0.9
  escalax=escala
  escalal=0.85
  ajuste=0.5
  escalat=0.5
  escalap=0.4
  listaFigurasGrafMtr=NULL
  contadorFiguras=1

  m.mat.global=NULL
  nombres=NULL
  i.min=NULL
  i.max=NULL
  x=NULL
  x.label=NULL
  x.name=NULL
  y.label=NULL
  m.global = pkg.env$m.global
  m.mat.global=m.global
  names.index = pkg.env$names.index
  i.min=1
  i.max=k.max-(k.min-1)
  names.metr = pkg.env$names.metr
  x=c(k.min:k.max)
  x.label="K values"
  x.name=as.character(c(k.min:k.max))
  y.label="Silhouette avg. width"
  #Pattern: QualityIndices__MetricX, ..., QualityIndices__MetricN
  figurename="QualityIndices_"

  i.min=k.min
  i.max=k.max
  margins <- par(mar=c(5,5,3,3))
  on.exit(par(margins))
  for (m.g in 1:length(names.metr)) {
    cur.k.width = m.mat.global[[m.g]][,1]
    cur.k.width = cur.k.width[!is.na(cur.k.width)]
    g.main=NULL
    leg.g=NULL
    ymin=NULL
    ymax=NULL
    xmin=NULL
    xmax=NULL
    xleg=NULL
    yleg=NULL
    c.max=NULL
    t.linea=NULL
    t.color=NULL
    ymarcas=NULL
    k.classes=NULL
    xmin=min(x)-0.25
    xmax=max(x)+0.25
    xleg=((xmax-xmin)*escalal)+3.2
    c.max=dim(m.mat.global[[m.g]])[2]
    ymin=min(cur.k.width)
    ymax=1
    ymarcas=round(seq(ymin,ymax,length.out=5),2)
    yleg=ymin+((ymax-ymin)/2)*seq(c.max,1,-1)/(2*c.max)
    t.linea=seq(1,c.max)
    t.color=rep("black",c.max)

    g.main=paste(" Qual. Indices of '", names.metr[m.g], "' for k in [",
                 k.min, ",", k.max,"]",sep="")

    # par(new=FALSE,bg="white",fg="black")
    y=NULL
    y.name=NULL
    y=cur.k.width
    y.name=names.index[1]
    leg.g[1] <- paste(y.name," avg. width",sep="")
    plot(x,y, type="l", xaxt="n", yaxt="n", xlab="", ylab="", main=g.main, xlim=c(xmin,xmax), ylim=c(ymin,ymax), lty=t.linea[1], col=t.color[1])
    par(new=TRUE)
    plot(x,y, type="o", xaxt="n", yaxt="n", xlab="", ylab="", main=g.main, xlim=c(xmin,xmax), ylim=c(ymin,ymax), lty=t.linea[1], col=t.color[1])
    par(new=TRUE)

    mtext(side=1, text=x.label,line=3)
    mtext(side=2, text=y.label,line=3)
    axis(side=1, at=x, labels=x.name, las=1, cex.axis=escalax)
    axis(side=2, at=ymarcas, labels=ymarcas, cex.axis=escalal)
    par(new=FALSE)
  }
}

runSilhouetteIMG <- function(data, k) {
  names.metr = pkg.env$names.metr
  datos.bruto = data
  estable = pkg.env$estable

  ancho=NULL
  alto=NULL
  ajuste=NULL
  figurename=NULL
  escala=NULL
  escalax=NULL
  escalal=NULL
  escalat=NULL
  escalap=NULL
  listaFigurasSil=NULL

  ancho=7
  alto=6
  escala=1     #new 0.6
  escalax=0.7  #new escala
  escalal=0.75  #new 0.8
  ajuste=0.5
  escalat=0.5
  escalap=0.4
  par(new=FALSE,bg="white",fg="black", cex=1, mex=.6)
  onto.matrix=NULL
  onto.matrix=matrix(data=NA, nrow=length(datos.bruto[,1]), ncol=(length(names.metr)+1))
  #onto.matrix[,1]=levels(datos.bruto[,1])
  onto.matrix[,1]=as.character(datos.bruto[,1])
  colnames(onto.matrix)=c("Datasets",paste(names.metr,sep="."))

  margenes=c(6,4,6,8)
  margins <- par(mar=margenes, cex=escala, mex=escalal)
  on.exit(par(margins))

  k.cl = k
  colores=c(2:(k.cl+1)) # 2 to k.cl+1, avoid number 1 since it's black and it's not pretty
  #Pattern: Silhouette_K_N_MetricX, ..., Silhouette_K_N_MetricN
  figurename="Silhouette_K_"
  contadorFiguras=1

  for (i.metr in 1:length(names.metr)) {
    datos=NULL
    i.datos=NULL
    name.file=NULL
    metric.onto=NULL
    metric.name=NULL
    part.onto=NULL
    sil.w=NULL

    sil.c=NULL
    x.leyenda=NULL
    t.leyenda=NULL
    xleyenda=NULL
    yleyenda=NULL
    x.leyenda=0.99

    metric.onto=datos.bruto[,i.metr+1]
    metric.name=names(datos.bruto)[i.metr+1]

    i.datos=i.metr
    if (estable[[i.datos]]$n.k==k.cl & estable[[i.datos]]$name.metric==metric.name) {
      part.onto=as.numeric(estable[[i.datos]]$kmk.dynamic.bs.or)
      onto.matrix[,(i.metr+1)]=part.onto

      sil.c=cluster.stats(dist(metric.onto), part.onto, noisecluster=TRUE)     #new

      sil.w=silhouette(part.onto, dist(metric.onto))
      estable[[i.datos]]$kmk.dynamic.bs.or.numeric=part.onto
      estable[[i.datos]]$sil.width=sil.w

      #DESDE AQUI new
      g.main=paste(metric.name,sep="")

      plot(sil.w, col=colores, main=g.main, border=NULL,
           mar=margenes, cex=escala, mex=escalal,
           cex.names = par("cex.axis"), do.n.k = TRUE, do.clus.stat = FALSE)
      t.leyenda=c(expression('j:  n'['j']), expression(' | ave'['i' %in% 'C'['j']]), expression('s'['i']))
      legend(x=x.leyenda,y=sil.c$n+1, legend=expression('j:  n'['j']), col="black",
             xjust=0, yjust=0, bty="n", xpd=TRUE, inset=c(-0.1,0), cex=escalax)
      legend(x=x.leyenda+0.03,y=sil.c$n+1, legend=expression(' | ave'['i' %in% 'C'['j']]), col="black",
             xjust=0, yjust=0, bty="n", xpd=TRUE, inset=c(-0.1,0), cex=escalax)
      legend(x=x.leyenda+0.1,y=sil.c$n+1, legend=expression(' s'['i']), col="black",
             xjust=0, yjust=0, bty="n", xpd=TRUE, inset=c(-0.1,0), cex=escalax)
      xleyenda=rep(x.leyenda,k.cl)
      yleyenda=(sil.c$cluster.size==1)*0.6*(sil.c$n-cumsum(sil.c$cluster.size))+
        (sil.c$n-cumsum(sil.c$cluster.size))+sil.c$cluster.size*3/k.cl+2
      leyenda=paste(names(sil.c$clus.avg.silwidths),rep(": ",k.cl),
                    sil.c$cluster.size,"|",round(sil.c$clus.avg.silwidths,digits=2),sec="")
      for (i.leyenda in 1:k.cl){
        legend(list(x=xleyenda[i.leyenda],y=yleyenda[i.leyenda]-0.7), legend=leyenda[i.leyenda], col="black",
               xjust=0, yjust=1, bty="n", xpd=TRUE, inset=c(-0.1,0), cex=escalal)
        }
    }
  }
}

runSilhouetteTable <- function(data, k) {
  names.metr = pkg.env$names.metr
  datos.bruto = data
  estable = pkg.env$estable
  k.cl = k
  ##
  #  Building table header
  ##
  silhouetteData <- list()
  silhouetteData$header <- list("Metric")
  for (i in 1:k.cl) {
    header = paste("Cluster_", i, "_SilScore", sep="")
    silhouetteData$header <- c(silhouetteData$header, header)
  }
  silhouetteData$header <- c(silhouetteData$header, "Avg_Silhouette_Width")
  for (i in 1:k.cl) {
    header = paste("Cluster_", i, "_Size", sep="")
    silhouetteData$header <- c(silhouetteData$header, header)
  }

  silhouetteData$header = unlist(silhouetteData$header, use.names=FALSE)

  ##
  #  //Building table header
  ##

  onto.matrix=NULL
  onto.matrix=matrix(data=NA, nrow=length(datos.bruto[,1]), ncol=(length(names.metr)+1))
  #onto.matrix[,1]=levels(datos.bruto[,1])
  onto.matrix[,1]=as.character(datos.bruto[,1])
  colnames(onto.matrix)=c("Datasets",paste(names.metr,sep="."))

  for (i.metr in 1:length(names.metr)) { # i.metr= n de metrica     i.metr=5
    datos=NULL
    i.datos=NULL
    name.file=NULL
    metric.onto=NULL
    metric.name=NULL
    part.onto=NULL
    sil.w=NULL
    sil.c=NULL
    x.leyenda=NULL
    xleyenda=NULL
    yleyenda=NULL
    metric.onto=datos.bruto[,i.metr+1]
    metric.name=names(datos.bruto)[i.metr+1]
    x.leyenda=0.99

    i.datos=i.metr
    if (estable[[i.datos]]$n.k==k.cl & estable[[i.datos]]$name.metric==metric.name) {
      part.onto=as.numeric(estable[[i.datos]]$kmk.dynamic.bs.or)
      onto.matrix[,(i.metr+1)]=part.onto
      sil.w=silhouette(part.onto, dist(metric.onto))
      sil.c=cluster.stats(dist(metric.onto), part.onto, noisecluster=TRUE)
      ## Building body rows
      silhouetteData$body[[i.metr]]=c(metric.name, sil.c$clus.avg.silwidths, mean(sil.w[,"sil_width"]), sil.c$cluster.size)
      silhouetteData$body[[i.metr]]=unlist(silhouetteData$body[[i.metr]], use.names=FALSE)
    }
  }  # end for i.metr

  silhouetteDataFrame = t(data.frame(silhouetteData$body))
  colnames(silhouetteDataFrame) = silhouetteData$header
  rownames(silhouetteDataFrame) <- NULL
  return(silhouetteDataFrame)
}

runSilhouetteTableRange <- function(data, k.min, k.max) {

  getHeader <- function(k) {
    silhouetteData$header <- list("Metric")
    for (i in 1:k) {
      header = paste("Cluster_", i, "_SilScore", sep="")
      silhouetteData$header <- c(silhouetteData$header, header)
    }
    silhouetteData$header <- c(silhouetteData$header, "Avg_Silhouette_Width")
    for (i in 1:k) {
      header = paste("Cluster_", i, "_Size", sep="")
      silhouetteData$header <- c(silhouetteData$header, header)
    }

    silhouetteData$header = unlist(silhouetteData$header, use.names=FALSE)
    return(silhouetteData$header)
  }

  names.metr = pkg.env$names.metr
  datos.bruto = data

  estable = pkg.env$estable
  k.min = k.min
  k.max = k.max

  onto.matrix=NULL
  onto.matrix=matrix(data=NA, nrow=length(datos.bruto[,1]), ncol=(length(names.metr)+1))
  onto.matrix[,1]=as.character(datos.bruto[,1])
  colnames(onto.matrix)=c("Datasets",paste(names.metr,sep="."))
  offset = 0
  estableLength = length(estable)
  names.metrLength = length(names.metr)

  silhouetteData <- list()
  silhouetteDataIndex <- vector(mode="integer", length=length(1:k.max))
  for (k in k.min:k.max) {
    header <- getHeader(k = k)
    silhouetteData[[k]] <- data.frame(matrix(ncol = length(header), nrow = names.metrLength))
    colnames(silhouetteData[[k]]) = header
    rownames(silhouetteData[[k]]) <- NULL
    silhouetteDataIndex[k] = 1
  }

  # estable object stores names.metr * length(k.min:k.max) entries
  for (i.metr in 1:estableLength) {
    cur.metr = as.integer(abs(i.metr-(names.metrLength*offset)))
    cur.data = estable[[i.metr]]
    cur.k = cur.data$n.k

    cur.row <- list(cur.data$name.metric)
    cur.row <- c(cur.row, cur.data$sil.c$clus.avg.silwidths)
    #cur.row <- c(cur.row, cur.data$sil.c$avg.silwidth)
    cur.row <- c(cur.row, mean(cur.data$sil.w[,"sil_width"]))
    cur.row <- c(cur.row, cur.data$sil.c$cluster.size)
    cur.row <- unlist(cur.row, use.names = FALSE)

    index = silhouetteDataIndex[cur.k]
    silhouetteData[[cur.k]] = insertRow(silhouetteData[[cur.k]], cur.row,index)
    silhouetteDataIndex[cur.k] = index + 1

    if (cur.metr == names.metrLength) { # Last metric
      offset = offset + 1
    }
  }

  ##
  # Data cleaning
  ##
  # Matrix inserts NA by default, remove them before returning the data
  for (k in k.min:k.max) {
    silhouetteData[[k]] <- na.omit(silhouetteData[[k]])
  }
  names(silhouetteData) <- paste("k_", 1:k.max, sep = "")
  silhouetteData[sapply(silhouetteData, is.null)] <- NULL
  return(silhouetteData)
}

checkKValue <- function(k) {
  if (k < 2 || k > 15) {
    error=paste("k value (",k,") is not in range [2,15]", sep="")
    stop(error)
  }
}