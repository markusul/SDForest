#' @export 
varImp <- function(object) {UseMethod("varImp")}

#' Extract variable importance of an SDTree
#' 
#' This function extracts the variable importance of an SDTree. 
#' The variable importance is calculated as the sum of the decrease in the loss 
#' function resulting from all splits that use this covariate.
#' @author Markus Ulmer
#' @param object an SDTree object
#' @return A named vector of variable importance
#' @seealso \code{\link{varImp.SDForest}} \code{\link{SDTree}}
#' @examples
#' data(iris)
#' tree <- SDTree(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width, iris)
#' varImp(tree)
#' @export
varImp.SDTree <- function(object){
  j_dec <- object$tree$Get(function(x)c(x$j, x$res_dloss), 
                           filterFun = function(x)!data.tree::isLeaf(x))
  var_importance <- rep(0, length(object$var_names))
  if(is.null(j_dec)){
    names(var_importance) <- object$var_names
    return(var_importance)
  }
  for(i in 1:ncol(j_dec)){
    var_importance[j_dec[1, i]] <- var_importance[j_dec[1, i]] + j_dec[2, i]
  }
  names(var_importance) <- object$var_names
  
  var_importance
}

#' Extract variable importance of an SDForest
#' 
#' This function extracts the variable importance of an SDForest.
#' The variable importance is calculated as the sum of the decrease in the loss function 
#' resulting from all splits that use a covariate for each tree. 
#' The mean of the variable importance of all trees results in the variable importance for the forest.
#' @author Markus Ulmer
#' @param object an SDForest object
#' @return A named vector of variable importance
#' @seealso \code{\link{varImp.SDTree}} \code{\link{SDForest}}
#' @aliases varImp
#' @examples
#' data(iris)
#' fit <- SDForest(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width, 
#'                  iris, nTree = 10)
#' varImp(fit)
#' @export
varImp.SDForest <- function(object){
  rowMeans(sapply(object$forest, varImp))
}