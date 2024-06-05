#' Print a SDTree
#' 
#' Print contents of the SDTree.
#' @author Markus Ulmer
#' @param object Fitted object of class \code{SDTree}.
#' @seealso \code{\link{SDTree}}
#' @export
print.SDTree <- function(object){
  # print function for the spectral deconfounded tree
  print(object$tree, 'value', 's', 'j', 'label', 'decision', 'n_samples')
}

#' Print SDForest
#' 
#' Print contents of the SDForest.
#' @author Markus Ulmer
#' @param object Fitted object of class \code{SDForest}.
#' @seealso \code{\link{SDForest}}
#' @export
print.SDForest <- function(object){
  cat("SDForest result\n\n")
  cat("Number of trees: ", length(object$forest), "\n")
  cat("Number of covariates: ", length(object$var_names), "\n")
  if(!is.null(object$oob_loss)){
    cat("OOB loss: ", round(object$oob_loss, 2), "\n")
    cat("OOB spectral loss: ", round(object$oob_SDloss, 2), "\n")
  }
}