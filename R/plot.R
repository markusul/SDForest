#' Plot SDTree
#' 
#' Plot the SDTree.
#' @author Markus Ulmer
#' @param x Fitted object of class \code{SDTree}.
#' @param ... Further arguments passed to or from other methods.
#' @seealso \code{\link{SDTree}}
#' @export
plot.SDTree <- function(x, ...){
  data.tree::SetEdgeStyle(x$tree, label = function(e) {e$decision})
  data.tree::SetNodeStyle(x$tree, label = function(n) {n$label})
  plot(x$tree)
}
