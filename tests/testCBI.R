library(evaluomeR)


evaluomeRSupportedCBI()

<<<<<<< HEAD
dataFrame <- stability(data=ontMetrics, cbi="kmeans", k=3, all_metrics=FALSE, bs=100, L1=10)
assay(dataFrame)

dataFrame <- stabilityRange(data=ontMetrics, cbi="rskc", k.range=c(3,4), all_metrics=TRUE, bs=100, L1=10)
assay(dataFrame)

dataFrame <- stabilitySet(data=ontMetrics, k.set=c(3,4), bs=100, cbi="rskc", all_metrics=TRUE, L1=10)
assay(dataFrame)

dataFrame <- quality(data=ontMetrics, cbi="rskc", k=3, all_metrics=TRUE, L1=10)
assay(dataFrame)

dataFrame <- qualityRange(data=ontMetrics, cbi="rskc", k.range=c(3,4), all_metrics=TRUE, L1=10)
assay(dataFrame$k_3)

dataFrame <- qualitySet(data=ontMetrics, cbi="rskc", k.set=c(3,5), all_metrics=TRUE, L1=10)
assay(dataFrame$k_3)


# RSKC will not work with a dataframe of 1 column

sim <-
  function(mu,f){
    D<-matrix(rnorm(60*f),60,f)
    D[1:20,1:50]<-D[1:20,1:50]+mu
    D[21:40,1:50]<-D[21:40,1:50]-mu
    return(D)
  }
sim
d0<-sim(1,500)# generate a dataset
true<-rep(1:3,each=20) # vector of true cluster labels
d<-d0
ncl<-3
for ( i in 1 : 10){
  d[sample(1:60,1),sample(1:500,1)]<-rnorm(1,mean=0,sd=15)
}

# The generated dataset looks like this...
pairs(
  d[,c(1,2,3,200)],col=true,
  labels=c("clustering feature 1",
           "clustering feature 2","clustering feature 3",
           "noise feature1"),
  main="The sampling distribution of 60 cases colored by true cluster labels",
  lower.panel=NULL)

d

# RSKC works when more than 2 columns are provided

r3<-RSKC(d[,1:5],ncl,alpha=10/60,L1=6,nstart=200)
=======

dataFrame <- stability(data=ontMetrics, cbi="kmeans", k=3, all_metrics=FALSE, bs=100)
assay(dataFrame)

dataFrame <- stabilityRange(data=ontMetrics, cbi="rskc", k.range=c(3,4), all_metrics=TRUE, bs=100, L1=2)
assay(dataFrame)

dataFrame <- stabilitySet(data=ontMetrics, k.set=c(3,4), bs=100, cbi="rskc", all_metrics=TRUE, L1=2)
assay(dataFrame)

dataFrame <- quality(data=ontMetrics, cbi="rskc", k=3, all_metrics=TRUE, L1=2)
assay(dataFrame)

dataFrame <- qualityRange(data=ontMetrics, cbi="rskc", k.range=c(3,4), all_metrics=TRUE, L1=2)
assay(dataFrame$k_3)

dataFrame <- qualitySet(data=ontMetrics, cbi="rskc", k.set=c(3,5), all_metrics=TRUE, L1=2)
assay(dataFrame$k_3)

>>>>>>> devel
