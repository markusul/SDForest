#' @export
copy <- function(object, ...) UseMethod('copy')

#' Copy a tree
#' 
#' Returns a copy of the tree object. 
#' Might be useful if you want to keep the original tree in comparison to the pruned tree.
#' @author Markus Ulmer
#' @param object A SDTree object
#' @return A copy of the SDTree object
#' @seealso \code{\link{prune}}
#' @aliases copy
#' @export
copy.SDTree <- function(object){
  new_tree <- data.tree::Clone(object$tree)
  new_object <- object
  new_object$tree <- new_tree
  return(new_object)
}

#' Copy a forest
#' 
#' Returns a copy of the forest object.
#' Might be useful if you want to keep the original forest in comparison to the pruned forest.
#' @author Markus Ulmer
#' @param object A SDForest object
#' @return A copy of the SDForest object
#' @seealso \code{\link{prune}}
#' @aliases copy
#' @export
copy.SDForest <- function(object){
  new_forest <- lapply(object$forest, function(tree){copy(tree)})
  new_object <- object
  new_object$forest <- new_forest
  return(new_object)
}