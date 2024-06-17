#' @title Dataframe getter for \code{qualityRange} function.
#' @name getDataQualityRange
#' @aliases getDataQualityRange
#' @description
#' This method is a wrapper to retrieve a specific \code{\link{SummarizedExperiment}} given a \code{k} value from
#' the object returned by \code{\link{qualityRange}} function.
#'
#' @param data The object returned by \code{\link{qualityRange}} function.
#' @param k The desired \code{k} cluster.
#'
#' @return The \code{\link{SummarizedExperiment}} that contains information about the selected \code{k} cluster.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' qualityRangeData <- qualityRange(ontMetrics, k.range=c(3,5), getImages = FALSE)
#' # Getting dataframe that contains information about k=5
#' k5Data = getDataQualityRange(qualityRangeData, 5)
#'
getDataQualityRange <- function(data, k) {
  dataNames = names(data)
  # From, e.g, k_4 to 4
  kValues = sapply(strsplit(dataNames, split='_', fixed=TRUE), function(x) (x[2]))
  kValues = as.integer(kValues)
  kValues.length = length(kValues)
  if (k >= kValues[1] && k <= kValues[kValues.length]) {
    column = paste("k_", k, sep="")

    return(data[[column]])
  } else {
    error=paste("Selected k (",k,") is not in the range of k.range ["
                , kValues[1], ",", kValues[kValues.length], "]", sep="")
    stop(error)
  }
}

###############################
## Package local environment/variables ##
###############################

pkg.env <- new.env()
pkg.env$m.stab.global = NULL
pkg.env$m.stab.global.csv = NULL
pkg.env$m.global = NULL
pkg.env$e.global = NULL
pkg.env$estable = NULL
pkg.env$names.index = NULL
pkg.env$names.metr = NULL
pkg.env$seed = 13606
pkg.env$cbi = c("kmeans", "clara", "clara_pam", "hclust", "pamk", "pamk_pam", "rskc") # Supported CBIs


#####################
## Private methods ##
#####################

# Source: http://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
# Supresses the output messages of a function 'x'
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

isString <- function(str) {
  if (is.character(str) & length(str) > 0) {
    return(TRUE)
  }
  error=paste("Input parameter '", str, "' is not a string", sep="")
  stop(error)
}

## Inserts a newrow to existingDF at position r
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

## Checks if k is in range [2,15]
checkKValue <- function(k) {
  if (k < 2 || k > 15) {
    stop("k value is not in range [2,15]")
  }
}

checkDirectory <- function(path) {
  if (is.null(path) || !dir.exists(path)) {
    error <- paste("Invalid path '", path, "'. Directory does not exists or it is not a directory", sep="")
    stop(error)
  }
  if (!endsWith(path, "/")) {
    path <- paste(path, "/", sep="")
  }
  return(path)
}

sample.range <- function(x) {diff(range(x))}
sample.min <- function(x) {min(x)}

# data: One dataframe, thus one assay
createSE <- function(data) {
  se = SummarizedExperiment(assays = data.matrix(data),
          colData = MultiAssayExperiment::DataFrame(metrics = colnames(data)))
  return(se)
}

# data: A list of SummarizedExperiment objects (ExperimentList)
createSEList <- function(data) {
  if (!is.list(data)) {
    stop("Input variable is not a list")
  }
  if (length(data) == 0) {
    stop("Input variable is an empty list")
  }
  length = length(names(data))
  seList <- list()
  for (i in 1:length) {
    cur.data <- data[[i]]
    dataMatrix <- as.matrix(cur.data)
    #cat("DataMatrix '", names(data)[i], "'\n")
    #print(dataMatrix)
    if (is.na(dataMatrix[1, "Metric"])) { # Metrics are NA? At least the first one
      dataMatrix[,1] <- cur.data$Metric
    }

    se <- createSE(dataMatrix)
    seList <- c(seList, se)
  }
  names(seList) <- names(data)
  expList <- ExperimentList(seList)
  return(expList)
}

## Returns a cluterboot interface (CBI) depending on the input string,
## with its specific function parameters in an list as in:
## list("method" = CBI Method, "args" = Lit of Specific arguments for this CBI)
helperGetCBI <- function(cbi=pkg.env$cbi, krange, ...) {
  cbi <- match.arg(cbi)
  inargs <- as.list(match.call(expand.dots = TRUE))
  switch(cbi,
         kmeans={
           args = list("krange" = krange)
           return(list("method" = kmeansCBI, "args" = args))
         },
         clara={
           args = list("k" = krange, "usepam"=FALSE)
           return(list("method" = claraCBI, "args" = args))
         },
         clara_pam={
           args = list("k" = krange, "usepam"=TRUE)
           return(list("method" = claraCBI, "args" = args))
         },
         hclust={
           args = list("k" = krange, "method"="ward.D2")
           return(list("method" = hclustCBI, "args" = args))
         },
         pamk={
           args = list("k" = krange, "usepam"=FALSE, criterion="asw")
           return(list("method" = pamkCBI, "args" = args))
         },
         pamk_pam={
           args = list("k" = krange, "usepam"=TRUE, criterion="asw")
           return(list("method" = pamkCBI, "args" = args))
         },
         rskc={
           args = list("k" = krange)
           if (!("L1" %in% names(inargs))) {
             print("No argument 'L1' provided. Computing best L1 boundry with 'sparcl::KMeansSparseCluster.permute', this might take a bit longer")
           }
           return(list("method" = rskcCBI, "args" = args))
         },
         {
           error=paste("Input CBI '", cbi, "' is not defined in the package", sep="")
           stop(error)
         }
  )
}
# Remove rows from a dataframe 'df' where for a certain column the row provides 'NA' value
removeNAValues <- function(df, verbose=TRUE) {
  rowNames = df[[1]]
  affected = which(rowSums(is.na(df)) > 0)
  if (verbose && length(affected) > 0) {
    message("Warning: There are rows with NA values. I will remove these...")
  }
  for (row_i in affected) {
    if (verbose) {
      message("Row '", rowNames[row_i], "' was removed. NA values found in columns:")
    }
    for (column in colnames(df)) {
      naDetected = is.na(df[row_i, column])
      if (verbose && naDetected) {
        message("- '", column,"' ")
      }
    }
  }
  df <- na.omit(df)
  if (verbose) {
    message("")
  }
  return (df)
}

# Basic stats of a dataframe
dfStats <- function(df) {
  message("Data loaded.\n",
          "Number of rows: ", nrow(df), "\n",
          "Number of columns: ", ncol(df), "\n\n"
  )
}

#' @title Get supported CBIs in evaluomeR.
#' @aliases evaluomeRSupportedCBI
#' @description
#' A vector of supported CBIs available in evaluomeR.
#'
#' @return A String vector.
#'
#' @examples
#' supportedCBIs <- evaluomeRSupportedCBI
#'
evaluomeRSupportedCBI <- function() {
  return(pkg.env$cbi)
}
