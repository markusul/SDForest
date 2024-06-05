#' Plot SDTree
#' 
#' Plot the SDTree.
#' @author Markus Ulmer
#' @param object Fitted object of class \code{SDTree}.
#' @seealso \code{\link{SDTree}}
plot.SDTree <- function(object){
  data.tree::SetEdgeStyle(object$tree, label = function(x) {x$decision})
  data.tree::SetNodeStyle(object$tree, label = function(x) {x$label})
  plot(object$tree)
}