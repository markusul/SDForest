#' @export 
toList <- function(object, ...) UseMethod('toList')

#' @export 
fromList <- function(object, ...) UseMethod('fromList')

#' SDTree toList method
#' 
#' Converts a the tree in a SDTree object from a 
#' This makes it substantially easier to save the tree to disk.
#' @author Markus Ulmer
#' @param tree A SDTree object
#' @return A list representation of the tree
#' @seealso \code{\link{fromList}}
#' @aliases toList
#' @export
toList.SDTree <- function(tree){
  tree$tree <- as.list(tree$tree)
  return(tree)
}

#' SDTree fromList method
#' 
#' Converts a list to a SDTree object.
#' 
fromList.SDTree <- function(tree){
  tree$tree <- data.tree::as.Node(tree$tree)
  return(tree)
}

toList.SDForest <- function(forest){
  forest$forest <- lapply(forest$forest, toList)
  return(forest)
}

fromList.SDForest <- function(forest){
  forest$forest <- lapply(forest$forest, fromList)
  return(forest)
}