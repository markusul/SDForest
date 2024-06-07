#' Partial dependence
#' 
#' This function calculates the partial dependence of a model on a single variable.
#' For that predictions are made for all observations in the dataset while varying 
#' the value of the variable of interest. The overall partial effect is the average
#' of all predictions. \insertCite{Friedman2001GreedyMachine}{SDForest}
#' @importFrom Rdpack reprompt
#' @references
#'   \insertAllCited{}
#' @author Markus Ulmer
#' @param object A model object that has a predict method that takes newdata as argument 
#' and returns predictions.
#' @param j The variable for which the partial dependence should be calculated.
#' Either the column index of the variable in the dataset or the name of the variable.
#' @param X The dataset on which the partial dependence should be calculated.
#' Should contain the same variables as the dataset used to train the model.
#' If NULL, tries to extract the dataset from the model object.
#' @param mc.cores Number of cores to use for parallel computation.
#' @return An object of class \code{partDependence} containing
#' \item{preds_mean}{The average prediction for each value of the variable of interest.}
#' \item{x_seq}{The sequence of values for the variable of interest.}
#' \item{preds}{The predictions for each value of the variable of interest for each observation.}
#' \item{j}{The name of the variable of interest.}
#' \item{xj}{The values of the variable of interest in the dataset.}
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' y <- sign(x) * 3 + rnorm(100)
#' model <- SDTree(x = x, y = y, Q_type = 'no_deconfounding')
#' pd <- partDependence(model, 1, X = x)
#' plot(pd)
#' @seealso \code{\link{SDForest}}, \code{\link{SDTree}}
#' @export
partDependence <- function(object, j, X = NULL, mc.cores = 1){
  j_name <- j
  if(is.character(j)){
    j <- which(names(X) == j)
  }
  
  if(is.null(X)){
    X <- object$X
    if(is.null(X)) stop('X must be provided if it is not part of the object')
    
  }
  X <- data.frame(X)
  
  if(!is.numeric(j)) stop('j must be a numeric or character')
  if(j > ncol(X)) stop('j must be smaller than p')
  if(j < 1) stop('j must be larger than 0')
  if(any(is.na(X))) stop('X must not contain missing values')
  
  x_seq <- seq(min(X[, j]), max(X[, j]), length.out = 100)
  
  if(mc.cores > 1){
    preds <- parallel::mclapply(x_seq, function(x){
      X_new <- X
      X_new[, j] <- x
      pred <- predict(object, newdata = X_new)
      return(pred)
    }, mc.cores = mc.cores)
  }else{
    preds <- pbapply::pblapply(x_seq, function(x){
      X_new <- X
      X_new[, j] <- x
      pred <- predict(object, newdata = X_new)
      return(pred)
    })
  }
  preds <- do.call(rbind, preds)
  preds_mean <- rowMeans(preds)
  
  res <- list(preds_mean = preds_mean, x_seq = x_seq, preds = preds, j = j_name, xj = X[, j])
  class(res) <- 'partDependence'
  return(res)
}

#' Plot partial dependence
#' 
#' This function plots the partial dependence of a model on a single variable.
#' @author Markus Ulmer
#' @param x An object of class \code{partDependence} returned by \code{\link{partDependence}}.
#' @param n_examples Number of examples to plot in addition to the average prediction.
#' @param ... Further arguments passed to or from other methods.
#' @return A ggplot object.
#' @seealso \code{\link{partDependence}}
#' @export
plot.partDependence <- function(x, n_examples = 19, ...){
  ggdep <- ggplot2::ggplot() + ggplot2::theme_bw()
  preds <- x$preds
  x_seq <- x$x_seq
  
  
  sample_examples <- sample(1:ncol(preds), min(n_examples, ncol(preds)))

  for(i in sample_examples){
    pred_data <- data.frame(x = x_seq, y = preds[, i])
    ggdep <- ggdep + ggplot2::geom_line(data = pred_data, 
                                        ggplot2::aes(x = x, y = y), col = 'grey')
  }
  
  ggdep <- ggdep + ggplot2::geom_line(data = data.frame(x = x_seq, y = x$preds_mean), 
                                      ggplot2::aes(x = x, y = y), col = '#08cbba', 
                                      linewidth = 1.5)
  ggdep <- ggdep + ggplot2::geom_rug(data = data.frame(x = x$xj, 
                                                       y = min(preds[, sample_examples])), 
                                     ggplot2::aes(x = x, y = y), 
                                     sides = 'b', col = '#949494')
  ggdep <- ggdep + ggplot2::ylab('f(x)') + ggplot2::ggtitle('Partial dependence')
  if(is.character(x$j)){
    ggdep <- ggdep + ggplot2::xlab(x$j)
  }else{
    ggdep <- ggdep + ggplot2::xlab(paste('x', x$j, sep = ''))
  }
  
  ggdep
}
