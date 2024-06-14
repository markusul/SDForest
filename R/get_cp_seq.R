#' @export
get_cp_seq <- function(object, ...) UseMethod('get_cp_seq')

#' Get the sequence of complexity parameters of an SDTree
#' 
#' This function extracts the sequence of complexity parameters of an SDTree that 
#' result in changes of the tree structure if pruned. Only cp values that differ
#' in the first three digits after the decimal point are returned.
#' @author Markus Ulmer
#' @param object an SDTree object
#' @param ... Further arguments passed to or from other methods.
#' @return A sequence of complexity parameters
#' @seealso \code{\link{regPath}} \code{\link{stabilitySelection}}
#' @export
get_cp_seq.SDTree <- function(object, ...){
  cp_seq <- unique(object$tree$Get('cp_max'))
  cp_seq[cp_seq > 1] <- 1
  cp_seq <- unique(ceiling(cp_seq * 1000)/1000)
  cp_seq <- c(0, cp_seq)
  return(cp_seq)
}

#' Get the sequence of complexity parameters of an SDForest
#' 
#' This function extracts the sequence of complexity parameters of an SDForest that
#' result in changes of the SDForest if pruned. Only cp values that differ
#' in the first three digits after the decimal point are returned.
#' @author Markus Ulmer
#' @param object an SDForest object
#' @param ... Further arguments passed to or from other methods.
#' @return A sequence of complexity parameters
#' @seealso \code{\link{regPath}} \code{\link{stabilitySelection}} 
#' \code{\link{get_cp_seq.SDTree}}
#' @aliases get_cp_seq
#' @export
get_cp_seq.SDForest <- function(object, ...){
  cp_seq <- unique((unlist(lapply(object$forest, function(x) x$tree$Get('cp_max')))))
  cp_seq[cp_seq > 1] <- 1
  cp_seq <- unique(ceiling(cp_seq * 1000)/1000)
  cp_seq <- c(0, cp_seq)
  return(cp_seq)
}