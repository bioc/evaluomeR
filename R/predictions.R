

# Funcion que calcula las componentes principales a partir de la matriz de las metricas de OQuaRE pasado por parametro. Imprime un
# vector que indica el Porcentaje de Varianza Explicado de cada una de las componentes y devuelve una matriz que representa
# las nuevas variables.
metrics_pca <- function(data){
  # Calculamos las componentes principales
  pca_result = prcomp(data, scale = TRUE)
  pca_result$rotation
  pca_result$x = - pca_result$x
  pca_result$rotation = - pca_result$rotation
  VE = pca_result$sdev^2
  PVE = VE / sum(VE) # Porcentaje de varianza explicado

  count = 1
  PVE_accumulated = PVE[count]
  aux = PVE[1] + PVE[2]
  while(count < length(data) & aux < 0.9){
    count = count + 1
    PVE_accumulated = PVE_accumulated + PVE[count]
    aux = aux + PVE[count+1]
  }
  message(paste0(count, " PC used, ", round(PVE_accumulated,2)*100, "% explained"))
  return (data.frame(pca_result$x[,1:count]))
}

metrics_randomforest <- function(data) {
  variables = randomForest(x=data)
  variables_importancia = variables$importance
  variables_importancia = as.data.frame(variables_importancia[order(variables_importancia, decreasing = TRUE),])
  colnames(variables_importancia) = colnames(variables$importance)

  most_important_metrics = head(rownames(variables_importancia), 7)
  return (data[, most_important_metrics])
}


# Funcion que devuelve el modelo que mejor se ajusta a los datos segun el criterio pasado por parametro. Para ello, define una
# serie de modelos, donde la diferencia entre cada uno de ellos sera el numero de componentes que contendra la mixtura, y
# finalmente se calcula cual es el mejor modelo.
# Parametros de entrada:
# - data: matriz de datos con la que se va a trabajar, en este caso va a ser la matriz de las componentes principales que
#         vamos a utilizar.
# - seed: semilla de inicializacion
# - k.range: rango de componentes para el cual vamos a definir un modelo.
# - nrep: numero de inicializaciones aleatorias utilizadas a la hora de ajustar los modelos.
# - criterion: criterio utilizado para seleccionar el mejor modelo. Admite los valores "BIC" y "AIC".

flemixModel = function(data, seed, k.range, nrep, criterion){
  old.seed = .Random.seed
  on.exit( { .Random.seed = old.seed } )

  model = list()
  nComponents = length(data)
  for (i in 1:nComponents){
    model[[i]] = FLXMRglm(data[,i]~.)
  }

  if (!is.null(seed)) {
    set.seed(as.integer(seed)) #seed
  } else {
    # Default seed
    set.seed(pkg.env$seed)
  }

  mixtures = stepFlexmix(~ 1, data = data, k = k.range, nrep = nrep, model = model)
  mixture_best = getModel(mixtures, criterion)
  return(mixture_best)
}

# Funcion que devuelve la prediccion del modelo para cada ontologia. La prediccion del modelo es el valor que creemos que
# nuestro modelo asignara a cada ontologia, y depende de la importancia que cada componente le asigna a cada ontologia y del
# valor que cada componente asigna a cada ontologia.
# Parametros de entrada:
# - data: Matriz de datos de entrada, en este caso son las componentes principales
# - model: Modelo que vamos a utilizar.
# Salida: matriz de tamano numero de ontolgias x numero de componentes del modelo, en donde para cada ontologia
#         se indica el valor que cada componente (distribucion) le asigna.

getPrediction = function(data, model){
  # Crear el dataframe que va a almacenar las predicciones
  M = ncol(data)
  N = nrow(data)
  numberComponents = model@k
  aux = numberComponents+1
  predictions = data.frame(matrix(nrow = N, ncol = aux))
  rownames(predictions) = rownames(data)
  col_names = c("Description")
  col_names = c(col_names, names(model@components))
  colnames(predictions) = col_names
  predictions[,1] = rownames(predictions)

  # Obtener los pesos que cada componente principal tiene para cada componente de la mixtura
  coef = flexmix::predict(model)
  coef_summary = list()
  for (i in 1:length(coef)){
    coef_summary[[i]] = coef[[i]][1,]
  }

  # Obtener las predicciones de cada ontologia
  for(i in 1:N){
    for (j in 1:numberComponents){
      aux = j+1
      predictions[i,aux] = sum(coef_summary[[j]]*data[i,])
    }

  }
  return (predictions)
}


# Funcion que calcula la metrica global sobre cada una de las ontolgias. Esta vendra definida por el valor (prediccion)
# que cada componente asigna a cada ontologia y el peso que esa componente tiene en el modelo.
# Parametros de entrada:
# - input:data: Matriz de datos de entrada, en este pueden ser las componentes principales o las metricas de OQuaRE.
# - PCA: valor booleano que indica si los datos de entrada son las componentes principales (PCA = TRUE) o si
#       son las metricas de OQuaRE y se deben calcuar las CP (PCA = FALSE).
# - seed: semilla de inicializacion.
# - k_range: rango de componentes para el cual vamos a definir un modelo.
# - nrep: numero de inicializaciones aleatorias utilizadas a la hora de ajustar los modelos.
# - criterion: criterio utilizado para seleccionar el mejor modelo. Admite los valores "BIC" y "AIC".

#' @title Global metric score defined by a prediction.
#' @name globalMetric
#' @aliases globalMetric
#' @description
#' This analysis calculates a global metric score based upon a prediction model
#' computed with \code{\link{flexmix}} package.
#'
#' @inheritParams stability
#' @param k.range Concatenation of two positive integers.
#' The first value \code{k.range[1]} is considered as the lower bound of the range,
#' whilst the second one, \code{k.range[2]}, as the higher. Both values must be
#' contained in [2,15] range.
#' @param nrep Positive integer. Number of random initializations used in adjusting the model.
#' @param PCA Boolean. If true, a PCA is performed on the input dataframe before computing
#' the predictions.
#' @param criterion String. Critirion applied in order to select the best model. Possible values:
#' "BIC" or "AIC".
#'
#' @return A dataframe containing the global metric score for each metric.
#'
#' @examples
#' # Using example data from our package
#' data("rnaMetrics")
#' globalMetric(rnaMetrics, k.range = c(2,3), nrep=10, criterion="AIC", PCA=TRUE)
#'
#'
#'
globalMetric = function(data, k.range=c(2,15), nrep=10, criterion=c("BIC", "AIC"), PCA=FALSE, seed=NULL) {

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
  if (nrep <= 0) {
    stop("Please 'nrep' parameter should be a positive integer > 1.")
  }

  data = as.data.frame(SummarizedExperiment::assay(data))
  rownames(data) = data[, 1]
  data = data[, 2:ncol(data)]

  # Si los datos de entrada no son las componentes principales, las calculamos
  if (PCA == TRUE){
    data = metrics_pca(data)
  } else {
    data = metrics_randomforest(data)
  }
  print(data)

  criterion =  match.arg(criterion)
  # Calcular el modelo
  mixture_best = flemixModel(data, seed, k.range, nrep, criterion)
  nComponents = mixture_best@k

  # Obtener las predicciones
  predictions = getPrediction(data, mixture_best)
  # Calcular la metrica global
  n = nrow(data)
  global_metric = data.frame(matrix(nrow = n, ncol = 2))
  global_metric[,1] = rownames(predictions)
  rownames(global_metric) = global_metric[,1]
  weigths = prior(mixture_best)
  for(i in 1:n){
    for (j in 1:nComponents){
      global_metric[i,2] = sum(weigths*as.numeric(predictions[i,-1]))
    }
  }
  colnames(global_metric) = c("Description", "Prediction")
  return (global_metric)
}
