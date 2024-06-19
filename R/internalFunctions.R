######################################################
#function
#   jaccard(x, y)
#     for two binary variables (values 0,1)

jaccard = function (x, y) {
  c11 = sum(x == 1 & y == 1)
  c10 = sum(x == 1 & y == 0)
  c01 = sum(x == 0 & y == 1)
  return (c11 / (c11 + c10 + c01))
}

######################################################
#function
#   jaccard.cluster(clust1, clust2)
#     for two clusterings from a data set

jaccard.cluster <- function(clust1,clust2) {
  cat <- sort(unique(c(clust1,clust2)))
  aux <- rep(0,length(cat))
  if (length(clust1)==length(clust2) & all(sort(unique(clust1))==sort(unique(clust2)))) {
    for (i.aux in cat) {
      clust.x <- clust1==i.aux
      for (j.aux in cat) {
        clust.y <- clust2==j.aux
        aux[i.aux] <- max(aux[i.aux],jaccard(clust.x,clust.y))
      }
    }
  }
  return(aux)
}

######################################################
#function
#   boot.cluster(data, nk=5, B=100, seed=NULL)
#     stability cluster based on Jaccard by clustering bootstrap

boot.cluster <- function(data, nk=5, B=10, seed=NULL, prnt=FALSE) {
  old.seed <- .Random.seed
  on.exit( { .Random.seed <- old.seed } )
  if (!is.null(seed)) {
    #http://stackoverflow.com/questions/14324096/setting-seed-locally-not-globally-in-r?rq=1
    set.seed(as.integer(seed)) #seed
  } else {
    # Default seed
    set.seed(pkg.env$seed)
  }
  data <- as.matrix(data)
  n.data <- nrow(data)
  cluster1 <- kmeans(x=data, centers=nk, iter.max=100)
  cluster1$partition <- cluster1$cluster
  cluster1$nk <- nk
  bs.jaccard <- matrix(0,nrow=nk, ncol=B)
  for (i in 1:B){  # for B botstrap replications
    if (prnt) cat("bootstrap replicate ", i,"\n")  # Show the bootstrap replication index
    bscheck <- TRUE
    while (bscheck) {
      bs.samp <- sample(n.data,n.data,replace=TRUE)
      bs.samp <- unique(bs.samp)
      bs.data <- data[bs.samp,]
      if (nk <= length(levels(as.factor(bs.data))))
        bscheck <- FALSE
    }
    if (length(bs.data) <= nk) { # Not enough data for kmeans
      message("\tWarning: Could not process data for k = ", nk)
      return (NULL)
    }

    bs.cluster <- kmeans(bs.data, centers=nk, iter.max=100)
    bs.cluster$partition <- bs.cluster$cluster
    bs.cluster$nk <- nk
    bs.jaccard[,i] <- jaccard.cluster(cluster1$partition[bs.samp], bs.cluster$partition)
  } # end for i in 1:B
  bs.jaccard.mean <- apply(bs.jaccard, 1, mean, na.rm=TRUE) # 1=by rows
  return(list(result=cluster1, partition=cluster1$partition, nk=cluster1$nk,
              B=B, bsjaccard=bs.jaccard, means=bs.jaccard.mean))
}
######################################################

######################################################
#function
#   fAnova(clusterResult, k, num.elements)
#     Custom anova computation: F function.

fAnova <- function(clusterResult, k, n) {
  #cat("k: ", k, "\n")
  #cat("clusterResult$betweenss: ", clusterResult$betweenss, "\n")
  #cat("clusterResult$tot.withinss: ", clusterResult$tot.withinss, "\n")
  #cat("(k-1): ", (k-1), "\n")
  #cat("(n-k): ", (n-k), "\n")
  #cat("d1: ", clusterResult$betweenss/(k-1), "\n")
  #cat("d2: ", (clusterResult$tot.withinss)/(n-k), "\n")
  #cat("r: ", (clusterResult$betweenss/(k-1))/(clusterResult$tot.withinss)/(n-k), "\n")
  #("r: ", (clusterResult$betweenss/(k-1))/((clusterResult$tot.withinss)/(n-k)), "\n")

  return ((clusterResult$betweenss/(k-1))/((clusterResult$tot.withinss)/(n-k)))
}
######################################################

######################################################
#function
#   getMeasureValue(km5, measureName)
#     Returns a measure as in km5$measureName

getMeasureValue <- function(km5, measureName) {
  if (startsWith(measureName, "cluster_")) {
    return (toString(unlist(km5$csv[measureName], use.names = FALSE)))
  } else {
    stop("Unknown measure '", measureName, "'. Stopping.")
  }
}
######################################################

######################################################
#function
#   clusterbootWrapper(data, B, bootmethod="boot",
#                     clustermethod=kmeansCBI, krange, seed, ...)
#     Wrapper method for clusterboot functionality.

clusterbootWrapper <- function(data, B, bootmethod="boot",
                               cbi, krange, seed, ...) {
  cbiHelperResult = helperGetCBI(cbi, krange, ...)

  #cat("Using: ", cbi, "\n")
  #cat("Type: ", typeof(cbiHelperResult[["method"]]), "\n")
  #print(cbiHelperResult)

  mandatoryArgs = list(
    "data"=data,
    "B"=B,
    "bootmethod"=bootmethod,
    "seed"=seed,
    "clustermethod"=cbiHelperResult[["method"]]
  )

  methodArgs = append(mandatoryArgs, cbiHelperResult[["args"]])
  # Append parameters that the user might have specified in the ellipsis
  methodArgs = append(methodArgs, list(...))

  return (
            #quiet(
              do.call(
                clusterboot,
                methodArgs
                #  )
            )
          )
}
######################################################

######################################################
#function
#   clusteringWrapper(data, cbi, krange, seed)
#     Wrapper method for clustering without bootstrap functionality.

clusteringWrapper <- function(data, cbi, krange, seed, ...) {
  cbiHelperResult = helperGetCBI(cbi, krange, ...)

  if(exists(".Random.seed")){
    # .Random.seed might not exist when launched as background job
    # so only store and reset if it exists
    old.seed <- .Random.seed
  }

  on.exit(
    {
      if(exists("old.seed")) {
        .Random.seed <<- old.seed
      }
    }
  )

  if (!is.null(seed)) set.seed(seed)


  #cat ("Using: ", cbi, "\n")

  mandatoryArgs = list(
    "data"=data
  )

  methodArgs = append(mandatoryArgs, cbiHelperResult[["args"]])
  # Append parameters that the user might have specified in the ellipsis
  methodArgs = append(methodArgs, list(...))

  # print(cbiHelperResult[["method"]])
  # print(methodArgs)

  return (
    quiet(
      do.call(
        cbiHelperResult[["method"]],
        methodArgs
      )
    )
  )
}
######################################################

######################################################
#function
#   checkIfCanCluster(data, k)
#       It prevents clusterboot method from getting stuck by checking if
#       there is enough data to perform a clustering.

checkIfCanCluster <- function(data, ...) {
  mc <- as.list(match.call(expand.dots = TRUE))
  k = NULL
  if ("krange" %in% names(mc)) {
    k = mc[["krange"]]
  } else if ("k" %in% names(mc) ) {
    k = mc[["k"]]
  } else {
    stop(paste0("Unexpected error. Could not retrieve 'k' value from dot-dot-dot argument. ",
                "Arguments were: ", names(mc)))
  }

  numUnique = length(unique(data))
  #print(paste0("Not unique:", length(data)))
  #print(paste0("Unique:", numUnique))
  #print(paste0("Division is:", numUnique/k))
  #print(paste0("Do I stop?: ", (numUnique/k) < 1))

  if ((numUnique/k) < 1) {
    stop(paste0("Not enough data to cluster '", numUnique, "' unique values in '", k, "' clusters"))
  }

}
######################################################
