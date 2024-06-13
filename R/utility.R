#' @importFrom Rdpack reprompt
#' @import GPUmatrix
#' @import DiagrammeR
#' @import data.tree
#' @importFrom stats lm.fit
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @importFrom stats predict
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats sd

#' @export
data.handler <- function(formula = NULL, data = NULL, x = NULL, y = NULL){
  if(is.null(formula)){
    if(is.null(x) | is.null(y)){
      stop("Error: Either data or x and y is required.")
    }else {
      if(is.vector(x)){
        x <- matrix(x, ncol = 1)
      }
      x <- apply(x, 2, function(x){if (is.character(x)) as.numeric(factor(x))
        else if(!is.numeric(x))  as.numeric(x)
        else x})
      if (!is.numeric(y)) stop("Error: y must be numeric. Only regression is supported at the moment.")
      if(any(is.na(x)) | any(is.na(y))){
        stop("Error: Missing values are not allowed.")
      }
      if(any(is.infinite(x)) | any(is.infinite(y))){
        stop("Error: Infinite values are not allowed.")
      }
      if(!is.numeric(y)){
        stop("Error: Only regression is suported at the moment. Y must be numeric.")
      }
      return(list(X = as.matrix(x), Y = as.numeric(y)))
    }
  }else {
    if(is.null(data)){
      stop("Error: data is required.")
    }else {
      Call <- match.call()
      indx <- match(c("formula", "data"), names(Call), nomatch = 0L)
      
      if (indx[1] == 0L) stop("a 'formula' argument is required")
      
      temp <- Call[c(1L, indx)]      # only keep the arguments we wanted
      temp[[1L]] <- quote(stats::model.frame) # change the function called
      m <- eval.parent(temp)
      
      Terms <- attr(m, "terms")
      if(any(attr(Terms, "order") > 1L)) stop("Trees cannot handle interaction terms")
      
      Y <- model.response(m)
      X <- model.matrix(attr(m, "terms"), m)[, -1L, drop = FALSE]
      
      if(any(is.infinite(X)) | any(is.infinite(Y))){
        stop("Error: Infinite values are not allowed.")
      }
      if(!is.numeric(Y)){
        stop("Error: Only regression is suported at the moment. Y must be numeric.")
      }
      list(X = X, Y = as.numeric(Y))
    }
  }
}

#' @export
predict_outsample <- function(tree, X){
  # predict for every observation in X f(x)
  # using the splitting rules from the tree
  if(is.null(dim(X))){
    return(traverse_tree(tree, X))
  }
  apply(X, 1, function(x)traverse_tree(tree, x))
}

#helper functions to label nodes for plotting
#' @export
splitt_names <- function(node, var_names = NULL){
  if(is.null(var_names)){
    node$label <- paste('X', node$j, ' <= ', round(node$s, 2), sep = '')
  }else{
    node$label <- paste(var_names[node$j], ' <= ', round(node$s, 2), sep = '')
  }
}
#' @export
leave_names <- function(node){
  new_name <- as.character(round(node$value, 1))
  if(new_name %in% node$Get('name', filterFun = data.tree::isLeaf)){
    new_name <- paste(new_name, '')
  }
  node$label <- new_name
}

# finds all the reasonable splitting points in a data matrix
#' @export
find_s <- function(X, max_candidates = 100){
  p <- ncol(X)
  if(p == 1){
    X <- matrix(X, ncol = 1)
  }
  n <- nrow(X)
  
  X_sort <- apply(X, 2, sort, method = 'quick')
  
  if(is.null(dim(X_sort))){
    X_sort <- matrix(X_sort, ncol = p)
  }
  
  # find middle points between observed x values
  s <- X_sort[-nrow(X_sort), ] + diff(X_sort)/2
  
  # for runtime reasons limit split candidates
  if(nrow(s) > max_candidates){
    s <- s[unique(floor(seq(1, dim(s)[1], length.out = max_candidates))), ]
  }
  
  if(is.null(dim(s))){
    matrix(s, ncol = p)
  }else{
    s
  }
}

#' @export
traverse_tree <- function(tree, x){
  # traverse the tree using the splitting rules and 
  # returns point estimate for f(x)
  if(tree$isLeaf){
    return(tree$value)
  }
  if(x[tree$j] <= tree$s){
    traverse_tree(tree$children[[1]], x)
  }else {
    traverse_tree(tree$children[[2]], x)
  }
}

#' @export
loss <- function(Y, f_X){
  as.numeric(sum((Y - f_X)^2) / length(Y))
}

#' @export
pruned_loss <- function(tree, X_val, Y_val, Q_val, t){
  # function to prune tree using the minimum loss decrease t
  # and return spectral loss on the validation set
  
  tree_t <- data.tree::Clone(tree)
  
  # prune tree
  data.tree::Prune(tree_t, function(x) x$dloss > t)
  
  # predict on test set
  f_X_hat_val <- predict_outsample(tree_t, X_val)
  
  # return spectral loss
  sum((Q_val %*% Y_val - Q_val %*% f_X_hat_val) ** 2) / length(Y_val)
}