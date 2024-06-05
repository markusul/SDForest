#' Predictions for the SDTree
#' 
#' Predicts the response for new data using a fitted SDTree.
#' @author Markus Ulmer
#' @param object Fitted object of class \code{SDTree}.
#' @param newdata New test data of class \code{data.frame} containing 
#' the covariates for which to predict the response.
#' @return A vector of predictions for the new data.
#' @examples
#' set.seed(1)
#' n <- 50
#' X <- matrix(rnorm(n * 20), nrow = n)
#' y <- sign(X[, 1]) * 3 + rnorm(n)
#' model <- SDTree(x = X, y = y, Q_type = 'no_deconfounding')
#' predict(model, newdata = data.frame(X))
#' @seealso \code{\link{SDTree}}
#' @export
predict.SDTree <- function(object, newdata){
  if(!is.data.frame(newdata)) stop('newdata must be a data.frame')
  if(!all(object$var_names %in% names(newdata))) stop('newdata must contain all covariates used for training')

  X <- newdata[, object$var_names]
  if(any(is.na(X))) stop('X must not contain missing values')
  return(predict_outsample(object$tree, X))
}

#' Predictions for the SDForest
#' 
#' Predicts the response for new data using a fitted SDForest.
#' @author Markus Ulmer
#' @param object Fitted object of class \code{SDForest}.
#' @param newdata New test data of class \code{data.frame} containing
#' the covariates for which to predict the response.
#' @return A vector of predictions for the new data.
#' @examples
#' set.seed(1)
#' n <- 50
#' X <- matrix(rnorm(n * 20), nrow = n)
#' y <- sign(X[, 1]) * 3 + rnorm(n)
#' model <- SDForest(x = X, y = y, Q_type = 'no_deconfounding', nTree = 10)
#' predict(model, newdata = data.frame(X))
#' @seealso \code{\link{SDForest}}
#' @export
predict.SDForest <- function(object, newdata){
  # predict function for the spectral deconfounded random forest
  # using the mean over all trees as the prediction
  # check data type
  if(!is.data.frame(newdata)) stop('newdata must be a data.frame')
  if(!all(object$var_names %in% names(newdata))) stop('newdata must contain all covariates used for training')

  X <- as.matrix(newdata[, object$var_names])
  if(any(is.na(X))) stop('X must not contain missing values')

  pred <- do.call(cbind, lapply(object$forest, function(x){predict_outsample(x$tree, X)}))
  return(rowMeans(pred))
}

#' Out-of-bag predictions for the SDForest
#' 
#' Predicts the response for the training data 
#' using only the trees in the SDForest 
#' that were not trained on the observation.
#' @author Markus Ulmer
#' @param object Fitted object of class \code{SDForest}.
#' @param X Covariates of the training data.
#' If \code{NULL}, the data saved in the object is used.
#' @return A vector of out-of-bag predictions for the training data.
#' @seealso \code{\link{SDForest}} \code{\link{prune.SDForest}} \code{\link{plotOOB}}
#' @export
predictOOB <- function(object, X = NULL){
  if(is.null(X)){
    X <- object$X
  }
  if(!is.matrix(X)) stop('X must be a matrix and either provided or in the object')

  n <- nrow(X)

  if(ncol(X) != length(object$var_names)){
    stop("The number of covariates has to correspond to the data used for training the forest.")
  }

  oob_ind <- object$oob_ind

  oob_predictions <- sapply(1:n, function(i){
    if(length(oob_ind[[i]]) == 0){
      return(NA)
    }
    xi <- X[i, ]
    predictions <- sapply(oob_ind[[i]], function(model){
      predict_outsample(object$forest[[model]]$tree, xi)
    })
    return(mean(predictions))
  })
  return(oob_predictions)
}