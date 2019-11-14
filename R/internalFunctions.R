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
  if (!is.null(seed)) {
    #http://stackoverflow.com/questions/14324096/setting-seed-locally-not-globally-in-r?rq=1
    old.seed <- .Random.seed
    on.exit( { .Random.seed <<- old.seed } )
    set.seed(as.integer(seed)) #seed
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
