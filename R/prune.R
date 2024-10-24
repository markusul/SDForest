#' @export
prune <- function(object, ...) UseMethod('prune')

#' Prune an SDTree
#' 
#' Removes all nodes that did not improve the loss by more than cp times the initial loss. 
#' Either by themselves or by one of their successors. Note that the tree is pruned in place.
#' If you intend to keep the original tree, make a copy of it before pruning.
#' @author Markus Ulmer
#' @param object an SDTree object
#' @param cp Complexity parameter, the higher the value the more nodes are pruned.
#' @param ... Further arguments passed to or from other methods.
#' @return A pruned SDTree object
#' @seealso \code{\link{copy}}
#' @export
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(50 * 20), nrow = 50)
#' Y <- rnorm(50)
#' tree <- SDTree(x = X, y = Y)
#' pruned_tree <- prune(copy(tree), 0.2)
#' tree
#' pruned_tree
prune.SDTree <- function(object, cp, ...){
  data.tree::Prune(object$tree, function(x) x$cp_max > cp)
  # new labels for leaves
  object$tree$Do(leave_names, filterFun = data.tree::isLeaf)
  # unknown new predictions
  object$predictions <- NULL
  # new variable importance
  object$var_importance <- varImp(object)
  
  object
}

#' Prune an SDForest
#' 
#' Prunes all trees in the forest and re-calculates the out-of-bag predictions and performance measures.
#' The training data is needed to calculate the out-of-bag statistics. Note that the forest is pruned in place.
#' If you intend to keep the original forest, make a copy of it before pruning.
#' @author Markus Ulmer
#' @param object an SDForest object
#' @param cp Complexity parameter, the higher the value the more nodes are pruned.
#' @param X The training data, if NULL the data from the forest object is used.
#' @param Y The training response variable, if NULL the data from the forest object is used.
#' @param Q The transformation matrix, if NULL the data from the forest object is used.
#' @param pred If TRUE the predictions are calculated, if FALSE only the out-of-bag statistics are calculated.
#' This can set to FALSE to save computation time if only the out-of-bag statistics are needed.
#' @param ... Further arguments passed to or from other methods.
#' @return A pruned SDForest object
#' @seealso \code{\link{copy}} \code{\link{prune.SDTree}} \code{\link{regPath}}
#' @aliases prune
#' @export
prune.SDForest <- function(object, cp, X = NULL, Y = NULL, Q = NULL, pred = TRUE, ...){
  pruned_forest <- lapply(object$forest, function(tree){prune(tree, cp)})
  object$forest <- pruned_forest

  if(is.null(X)) X <- object$X
  if(is.null(Y)) Y <- object$Y
  if(is.null(Q)) Q <- object$Q
  if(is.null(X) | is.null(Y)){
    stop('X and Y must either be provided or in the object')
  }

  n <- length(Y)

  if(is.null(Q)) {
    Q <- function(v) v
    warning('Q was not provided, using Identity matrix')
  }

  if(any(c(nrow(X), n) != length(object$predictions))){
    stop("The data has to correspond to the data used for training the forest.")
  }
  if(ncol(X) != length(object$var_names)){
    stop("The number of covariates has to correspond to the data used for training the forest.")
  }

  oob_predictions <- predictOOB(object, X)
  object$oob_predictions <- oob_predictions
  object$oob_SDloss <- loss(Q(Y), Q(oob_predictions))
  object$oob_loss <- loss(Y, oob_predictions)

  if(pred){
    # predict with all trees
    pred <- do.call(cbind, lapply(object$forest, function(x){matrix(predict_outsample(x$tree, X))}))
      
    # use mean over trees as final prediction
    f_X_hat <- rowMeans(pred)
    object$predictions <- f_X_hat
  }else {
    object$predictions <- NULL
  }
  # variable importance
  object$var_importance <- rowMeans(as.matrix(sapply(object$forest, function(x){matrix(x$var_importance)})))  

  object
}
