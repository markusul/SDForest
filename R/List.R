#' @export 
toList <- function(object, ...) UseMethod('toList')

#' @export 
fromList <- function(object, ...) UseMethod('fromList')

#' SDTree toList method
#' 
#' Converts the tree in an SDTree object from 
#' class \code{Node} \insertCite{Glur2023Data.tree:Structure}{SDForest} to class \code{list}.
#' This makes it substantially easier to save the tree to disk.
#' @author Markus Ulmer
#' @references
#'  \insertAllCited{}
#' @param object an SDTree object with the tree in Node format
#' @param ... Further arguments passed to or from other methods.
#' @return an SDTree object with the tree in list format
#' @seealso \code{\link{fromList}}
#' @export
toList.SDTree <- function(object, ...){
  object$tree <- as.list(object$tree)
  object
}

#' SDTree fromList method
#' 
#' Converts the tree in an SDTree object from
#' class \code{list} to class \code{Node} \insertCite{Glur2023Data.tree:Structure}{SDForest}.
#' @author Markus Ulmer
#' @references
#'  \insertAllCited{}
#' @param object an SDTree object with the tree in list format
#' @param ... Further arguments passed to or from other methods.
#' @return an SDTree object with the tree in Node format
#' @seealso \code{\link{toList}}
#' @export
fromList.SDTree <- function(object, ...){
  object$tree <- data.tree::as.Node(object$tree)
  object
}

#' SDForest toList method
#' 
#' Converts the trees in an SDForest object from
#' class \code{Node} \insertCite{Glur2023Data.tree:Structure}{SDForest} to class \code{list}.
#' This makes it substantially easier to save the forest to disk.
#' @author Markus Ulmer
#' @references
#'  \insertAllCited{}
#' @param object an SDForest object with the trees in Node format
#' @param ... Further arguments passed to or from other methods.
#' @return an SDForest object with the trees in list format
#' @seealso \code{\link{fromList}} \code{\link{toList.SDTree}}
#' @aliases toList
#' @export
toList.SDForest <- function(object, ...){
  object$forest <- lapply(object$forest, toList)
  object
}

#' SDForest fromList method
#' 
#' Converts the trees in an SDForest object from
#' class \code{list} to class \code{Node} \insertCite{Glur2023Data.tree:Structure}{SDForest}.
#' @author Markus Ulmer
#' @references
#'  \insertAllCited{}
#' @param object an SDForest object with the trees in list format
#' @param ... Further arguments passed to or from other methods.
#' @return an SDForest object with the trees in Node format
#' @seealso \code{\link{fromList}} \code{\link{fromList.SDTree}}
#' @aliases fromList
#' @export
fromList.SDForest <- function(object, ...){
  object$forest <- lapply(object$forest, fromList)
  object
}