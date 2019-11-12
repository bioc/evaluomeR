#' @title Calculation of Pearson correlation coefficient.
#' @name metricsCorrelations
#' @aliases metricsCorrelations
#' @description
#' Calculation of Pearson correlation coefficient between every pair of
#' metrics available in order to quantify their interrelationship degree.
#' The score is in the range [-1,1]. Perfect correlations: -1 (inverse),
#' and 1 (direct).
#'
#' @inheritParams stability
#' @param margins See \code{\link{par}}.
#'
#' @return The Pearson correlation matrix as an assay
#' in a \code{\link{SummarizedExperiment}} object.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' cor = metricsCorrelations(ontMetrics, getImages = TRUE, margins = c(1,0,5,11))
#'
metricsCorrelations <- function(data, margins=c(0,10,9,11), getImages=TRUE) {

  data <- as.data.frame(assay(data))

  MatCorr <- cor(data[,2:length(data)])

  if (getImages == TRUE) {
    runMetricsCorrelationIMG(data, MatCorr, margins)
  }
  se <- createSE(MatCorr)
  return(se)
}

runMetricsCorrelationIMG <- function(data, correlations, margins) {

  datos.bruto = data
  MatCorr = correlations
  names.metr=names(datos.bruto)[-c(1)]  #metric names

  ##########################################################
  ancho=5 #dimension de grafica  (pulgadas)
  alto=4 #dimension de grafica  (pulgadas)
  escala=0.6 #reescalamiento de texto
  escalax=escala #reescalamiento de ejes
  escalal=0.8 #reescalamiento de etiquetas de ejes
  ajuste=0.5 #ajuste 0=izq, 0.5=centrado(defecto), 1=derecha
  escalat=0.5
  escalap=0.4
  #margenes=c(1,0,5,11)
  margenes=margins
  margins <- par(mar=margenes)
  on.exit(par(margins))
  ## para graficos de correlaciones
  col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white",
                             "cyan", "#007FFF", "blue","#00007F"))
  #col1 <- colorRampPalette(c("black", "grey70"))
  col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                             "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
  col3 <- colorRampPalette(c("red", "white", "blue"))
  col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F",
                             "cyan", "#007FFF", "blue","#00007F"))
  wb <- c("white","black")
  ##########################################################

  #Pattern: Correlations_X_metrics
  figurename="Correlations_"

  ##########################################################

  ## circle + colorful number
  corrplot(MatCorr, is.corr = TRUE,
           method = "number", #"circle", "square", "ellipse", "number", "shade", "color", "pie"
           order = "alphabet",  #"original", "AOE", "FPC", "hclust", "alphabet"
           diag = FALSE, type = "lower", #"full", "lower", "upper"
           tl.pos = "n", #pos text labels: lt", "ld", "td", "d", "n"
           tl.cex = 0.8, tl.col = "black", tl.offset = 0.4, tl.srt = 90, #ratio, col, pos & rotation
           col = col1(200), outline=FALSE, #col of symbol & outline of symbol
           cl.pos = "n", #pos col labels: "r", "b", "n"
           number.cex = 0.7,  #ratio to write corr in table
           number.font = 1, #font to write corr in table
           number.digits = 2,#NULL, #number of digits to write corr in table
           win.asp = 1, #ratio of symbol in table
           mar=margenes, #margenes del grafico
           bg = "white", title = "", fg="black",
           addgrid.col = NULL, #col of grid: NULL, NA, "red", ...
           addCoef.col = NULL, #col of coefficients: NULL, "red", ...
           addCoefasPercent = FALSE, #coefficients in %
           add = FALSE)

  corrplot(MatCorr, is.corr = TRUE,
           method = "circle", #"circle", "square", "ellipse", "number", "shade", "color", "pie"
           order = "alphabet",  #"original", "AOE", "FPC", "hclust", "alphabet"
           diag = TRUE, type = "upper", #"full", "lower", "upper"
           tl.pos = "t", #pos text labels: lt", "ld", "td", "d", "n"
           tl.cex = 0.8, tl.col = "black", #tl.offset = 0.4, tl.srt = 90, #ratio, col, pos & rotation
           col = col1(200), outline=FALSE, #col of symbol & outline of symbol
           cl.pos = "r", #pos col labels: "r", "b", "n"
           cl.cex = 0.8, cl.ratio = 0.15, cl.align.text = "r", cl.offset = -0.25,
           win.asp = 1, #ratio of symbol in table
           bg = "white", title = "", fg="black",
           #cl.pos = "n", #pos col labels: "r", "b", "n"
           #mar=margenes, #margenes del grafico
           add = TRUE)
}
