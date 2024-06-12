#' @export 
regPath <- function(object, ...) UseMethod('regPath')

#' Calculate the regularization path of a SDTree
#'
#' This function calculates the variable importance of a SDTree
#' for different complexity parameters.
#' @author Markus Ulmer
#' @param object A SDTree object
#' @param cp_seq A sequence of complexity parameters.
#' If NULL, the sequence is calculated automatically using only relevant values.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class \code{paths} containing
#' \item{cp}{The sequence of complexity parameters.}
#' \item{varImp_path}{A \code{matrix} with the variable importance
#' for each complexity parameter.}
#' @seealso \code{\link{plot.paths}} \code{\link{prune}} \code{\link{get_cp_seq}} \code{\link{SDTree}}
#' @examples
#' set.seed(1)
#' n <- 10
#' X <- matrix(rnorm(n * 5), nrow = n)
#' y <- sign(X[, 1]) * 3 + sign(X[, 2]) + rnorm(n)
#' model <- SDTree(x = X, y = y, Q_type = 'no_deconfounding')
#' paths <- regPath(model)
#' plot(paths)
#' \dontrun{
#' plot(paths, plotly = TRUE)
#' }
#' @export
regPath.SDTree <- function(object, cp_seq = NULL, ...){
  object$tree <- data.tree::Clone(object$tree)
  
  if(is.null(cp_seq)) cp_seq <- get_cp_seq(object)
  cp_seq <- sort(cp_seq)

  res <- lapply(cp_seq, function(cp){
    pruned_object <- prune(object, cp)
    return(list(var_importance = pruned_object$var_importance))})

  varImp_path <- t(sapply(res, function(x)x$var_importance))
  colnames(varImp_path) <- object$var_names

  paths <- list(cp = cp_seq, varImp_path = varImp_path)
  class(paths) <- 'paths'
  
  paths
}

#' Calculate the regularization path of a SDForest
#' 
#' This function calculates the variable importance of a SDForest
#' and the out-of-bag performance for different complexity parameters.
#' @author Markus Ulmer
#' @param object A SDForest object
#' @param cp_seq A sequence of complexity parameters.
#' If NULL, the sequence is calculated automatically using only relevant values.
#' @param X The training data, if NULL the data from the forest object is used.
#' @param Y The training response variable, if NULL the data from the forest object is used.
#' @param Q The transformation matrix, if NULL the data from the forest object is used.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class \code{paths} containing
#' \item{cp}{The sequence of complexity parameters.}
#' \item{varImp_path}{A \code{matrix} with the variable importance
#' for each complexity parameter.}
#' \item{loss_path}{A \code{matrix} with the out-of-bag performance
#' for each complexity parameter.}
#' \item{cp_min}{The complexity parameter with the lowest out-of-bag performance.}
#' @seealso \code{\link{plot.paths}} \code{\link{plotOOB}} \code{\link{regPath.SDTree}} \code{\link{prune}} \code{\link{get_cp_seq}} \code{\link{SDForest}}
#' @aliases regPath
#' @examples
#' set.seed(1)
#' n <- 10
#' X <- matrix(rnorm(n * 5), nrow = n)
#' y <- sign(X[, 1]) * 3 + sign(X[, 2]) + rnorm(n)
#' model <- SDForest(x = X, y = y, Q_type = 'no_deconfounding')
#' paths <- regPath(model)
#' plotOOB(paths)
#' plot(paths)
#' \dontrun{
#' plot(paths, plotly = TRUE)
#' }
#' @export
regPath.SDForest <- function(object, cp_seq = NULL, X = NULL, Y = NULL, Q = NULL, ...){
  if(is.null(cp_seq)) cp_seq <- get_cp_seq(object)
  cp_seq <- sort(cp_seq)
  object$forest <- lapply(object$forest, function(tree){
    tree$tree <- data.tree::Clone(tree$tree)
    return(tree)
    })

  res <- pbapply::pblapply(cp_seq, function(cp){
    pruned_object <- prune(object, cp, X, Y, Q, pred = FALSE)
    return(list(var_importance = pruned_object$var_importance, 
                oob_SDloss = pruned_object$oob_SDloss, 
                oob_loss = pruned_object$oob_loss))})

  varImp_path <- t(sapply(res, function(x)x$var_importance))
  colnames(varImp_path) <- object$var_names

  loss_path <- t(sapply(res, function(x) c(x$oob_SDloss, x$oob_loss)))
  colnames(loss_path) <- c('oob SDE', 'oob MSE')
  paths <- list(cp = cp_seq, varImp_path = varImp_path, loss_path = loss_path,
                cp_min = cp_seq[which.min(loss_path[, 1])])
  class(paths) <- 'paths'
  
  paths
}

#' @export
stabilitySelection <- function(object, ...) UseMethod('stabilitySelection')

#' Calculate the stability selection of a SDForest
#' 
#' This function calculates the stability selection of a SDForest
#' \insertCite{Meinshausen2010StabilitySelectionb}{SDForest}.
#' Stability selection is calculated as the fraction of trees in the forest
#' that select a variable for a split at each complexity parameter.
#' @importFrom Rdpack reprompt
#' @references
#'   \insertAllCited{}
#' @author Markus Ulmer
#' @param object A SDForest object
#' @param cp_seq A sequence of complexity parameters.
#' If NULL, the sequence is calculated automatically using only relevant values.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class \code{paths} containing
#' \item{cp}{The sequence of complexity parameters.}
#' \item{varImp_path}{A \code{matrix} with the stability selection
#' for each complexity parameter.}
#' @seealso \code{\link{plot.paths}} \code{\link{regPath}} \code{\link{prune}} \code{\link{get_cp_seq}} \code{\link{SDForest}}
#' @aliases stabilitySelection
#' @examples
#' set.seed(1)
#' n <- 10
#' X <- matrix(rnorm(n * 5), nrow = n)
#' y <- sign(X[, 1]) * 3 + sign(X[, 2]) + rnorm(n)
#' model <- SDForest(x = X, y = y, Q_type = 'no_deconfounding')
#' paths <- stabilitySelection(model)
#' plot(paths)
#' \dontrun{
#' plot(paths, plotly = TRUE)
#' }
#' @export
stabilitySelection.SDForest <- function(object, cp_seq = NULL, ...){
  if(is.null(cp_seq)) cp_seq <- get_cp_seq(object)
  cp_seq <- sort(cp_seq)

  imp <- pbapply::pblapply(object$forest, function(x)regPath(x, cp_seq)$varImp_path > 0)

  imp <- lapply(imp, function(x)matrix(as.numeric(x), ncol = ncol(x)))
  imp <- Reduce('+', imp) / length(object$forest)
  colnames(imp) <- object$var_names
  paths <- list(cp = cp_seq, varImp_path = imp)
  class(paths) <- 'paths'
  
  paths
}

#' Visualize the paths of a SDTree or SDForest
#' 
#' This function visualizes the variable importance of a SDTree or SDForest
#' for different complexity parameters. Both the regularization path and
#' the stability selection path can be visualized.
#' @author Markus Ulmer
#' @param x A \code{paths} object
#' @param plotly If TRUE the plot is returned interactive using plotly. Might be slow for large data.
#' @param selection A vector of indices of the covariates to be plotted. 
#' Can be used to plot only a subset of the covariates in case of many covariates.
#' @param log_scale If TRUE the y-axis is on a log scale.
#' @param ... Further arguments passed to or from other methods.
#' @return A \code{ggplot} object with the variable importance for different regularization.
#' If the \code{path} object includes a cp_min value, a black dashed line is
#' added to indicate the out-of-bag optimal variable selection.
#' @seealso \code{\link{regPath}} \code{\link{stabilitySelection}}
#' @export
plot.paths <- function(x, plotly = FALSE, selection = NULL, log_scale = FALSE, ...){
  varImp_path <- x$varImp_path
  if(!is.null(selection)){
    varImp_path <- varImp_path[, selection]
  }

  imp_data <- data.frame(varImp_path, cp = x$cp)
  imp_data <- tidyr::gather(imp_data, key = 'covariate', value = 'importance', -cp)

  if(log_scale){
    imp_data$importance <- log(imp_data$importance + 1)
  }
  
  gg_path <- ggplot2::ggplot(imp_data, ggplot2::aes(x = cp, y = importance, 
                                                    col = covariate)) +
      ggplot2::geom_line() + 
      ggplot2::theme_bw() + 
      ggplot2::geom_rug(data = imp_data, ggplot2::aes(x = cp, y = importance), 
                        sides = 'b', col = '#949494')

  if(!is.null(x$cp_min)){
    gg_path <- gg_path + ggplot2::geom_vline(xintercept = x$cp_min, linetype = 'dashed')
  }

  if(plotly){
    return(plotly::ggplotly(gg_path))
  }else if(length(unique(imp_data$covariate)) > 20){
    gg_path + ggplot2::theme(legend.position = 'none')
  }else{
    gg_path
  }
}

#' Visualize the out-of-bag performance of a SDForest
#' 
#' This function visualizes the out-of-bag performance of a SDForest
#' for different complexity parameters. Can be used to choose the optimal
#' complexity parameter.
#' @author Markus Ulmer
#' @param object A paths object with loss_path \code{matrix} 
#' with the out-of-bag performance for each complexity parameter.
#' @return A ggplot object
#' @seealso \code{\link{regPath.SDForest}}
#' @export
plotOOB <- function(object){
    loss_data <- data.frame(object$loss_path, cp = object$cp)
    gg_sde <- ggplot2::ggplot(loss_data, ggplot2::aes(x = cp, y = oob.SDE)) +
        ggplot2::geom_line() + 
        ggplot2::theme_bw()

    gg_mse <- ggplot2::ggplot(loss_data, ggplot2::aes(x = cp, y = oob.MSE)) +
        ggplot2::geom_line() + 
        ggplot2::theme_bw()
    gridExtra::grid.arrange(gg_sde, gg_mse, ncol = 2)
}
