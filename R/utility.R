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
#' @importFrom utils lsf.str
#' @importFrom stats rbeta


data.handler <- function(formula = NULL, data = NULL, x = NULL, y = NULL){
  if(is.null(formula)){
    if(is.null(x) | is.null(y)){
      stop("Error: Either data or x and y is required.")
    }else {
      if(is.vector(x)) x <- as.matrix(x)
      if(!is.matrix(x)){
        stop("Error: x must be a matrix!")
      }
      if(!is.numeric(x)){
        stop("Error: x must contain numerical values! Use formula for categorical predictors.")
      }
      if (!is.numeric(y)) stop("Error: y must be numeric. Only regression is supported at the moment.")
      if(any(is.na(x)) | any(is.na(y))){
        stop("Error: Missing values are not allowed.")
      }
      if(any(is.infinite(x)) | any(is.infinite(y))){
        stop("Error: Infinite values are not allowed.")
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
      
      # ordinal covariates to numeric
      ord <- names(m)[sapply(m, is.ordered)]
      m[ord] <- lapply(m[ord], as.integer)
      
      Terms <- attr(m, "terms")
      if(any(attr(Terms, "order") > 1L)) stop("Trees cannot handle interaction terms")
      
      Y <- model.response(m)
      X <- model.matrix(attr(m, "terms"), m)[, -1L, drop = FALSE]
      
      if(any(is.infinite(X)) | any(is.infinite(Y))){
        stop("Error: Infinite values are not allowed.")
      }
      if(!is.null(Y) & !is.numeric(Y)){
        stop("Error: Only regression is suported at the moment. Y must be numeric.")
      }
      list(X = X, Y = Y)
    }
  }
}

predict_outsample <- function(tree, X){
  # predict for every observation in X f(x)
  # using the splitting rules from the tree
  if(is.null(dim(X))){
    return(traverse_tree(tree, X))
  }
  apply(X, 1, function(x)traverse_tree(tree, x))
}

#helper functions to label nodes for plotting

split_names <- function(node, var_names = NULL){
  if(is.null(var_names)){
    node$label <- paste('X', node$j, ' <= ', round(node$s, 2), sep = '')
  }else{
    node$label <- paste(var_names[node$j], ' <= ', round(node$s, 2), sep = '')
  }
}

leave_names <- function(node){
  new_name <- as.character(round(node$value, 1))
  if(new_name %in% node$Get('name', filterFun = data.tree::isLeaf)){
    new_name <- paste(new_name, '')
  }
  node$label <- new_name
}

# finds all the reasonable spliting points in a data matrix
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

loss <- function(Y, f_X){
  as.numeric(sum((Y - f_X)^2) / length(Y))
}

pruned_loss <- function(tree, X_val, Y_val, Q_val, t){
  # function to prune tree using the minimum loss decrease t
  # and return spectral loss on the validation set
  
  tree_t <- data.tree::Clone(tree)
  
  # prune tree
  data.tree::Prune(tree_t, function(x) x$dloss > t)
  
  # predict on test set
  f_X_hat_val <- predict_outsample(tree_t, X_val)
  
  # return spectral loss
  sum((Q_val(Y_val) - Q_val(f_X_hat_val)) ** 2) / length(Y_val)
}

# more efficient transformations
get_Qf <- function(X, type, trim_quantile = 0.5, q_hat = 0, gpu = FALSE, scaling = TRUE){
  if(gpu) ifelse(GPUmatrix::installTorch(), 
                 gpu_type <- 'torch', 
                 gpu_type <- 'tensorflow')

  if(type == 'no_deconfounding') {
    return(function(v) v)
  }
  X <- scale(X, center = TRUE, scale = scaling)
  
  svd_error <- function(X, q, f = 1, count = 1){
    tryCatch({
      svd(X * f, nv = 0, nu = q)
    }, error = function(e) {
      warning(paste(e, ':X multipied by number close to 1'))
      if(count > 5) stop('svd did not converge')
      return(svd_error(X, q, 1 + 0.0000000000000001 * 10 ^ count, count + 1))})
  }
  
  if(ncol(X) == 1){
    warning('only one covariate, no deconfounding possible')
    return(function(v) v)
  }
  
  modes <- c('trim' = 1, 'pca' = 2, 'no_deconfounding' = 3)
  if(!(type %in% names(modes))) stop(paste("type must be one of:", 
                                           paste(names(modes), collapse = ', ')))
  
  # number of observations
  n <- dim(X)[1]
  
  # needed number of singular values
  q <- q_hat
  if(type == 'trim'){
    q <- floor(quantile(1:min(dim(X)), 1-trim_quantile))
  }
  
  # calculate deconfounding matrix
  sv <- svd_error(X, q)
  Uq <- sv$u[, 1:q]

  if(gpu) Uq <- gpu.matrix(Uq, type = gpu_type)

  switch(modes[type], 
         {#trim
           D_tilde <- sv$d[1:q]
           D_tilde <- D_tilde[q] / D_tilde
           
           #tau <- quantile(sv$d, trim_quantile)
           #D_tilde <- unlist(lapply(sv$d, FUN = function(x)min(x, tau))) / sv$d
           #D_tilde[is.na(D_tilde)] <- 1
           #q <- sum(D_tilde != 1)
           #D_tilde <- D_tilde[1:q]
           },
         {# pca
           if(q_hat <= 0) 
             stop("the assumed confounding dimension must be larger than zero, increase q_hat")
           D_tilde <- rep(0, q_hat)
           #q <- q_hat
           }
         )
  
  
  Qf <- function(v){
    UqV <- crossprod(Uq, v)
    v + Uq %*% (UqV * (D_tilde - 1))
  }
  return(Qf)
}

get_Wf <- function(A, gamma, intercept = FALSE, gpu = FALSE){
  if(intercept) A <- cbind(1, A)
  if(ncol(A) > nrow(A)) stop('A must have full rank!')
  if(gamma < 0) stop('gamma must be non-negative')
  
  if(gpu) ifelse(GPUmatrix::installTorch(), 
                 gpu_type <- 'torch', 
                 gpu_type <- 'tensorflow')  
  if(gpu) A <- gpu.matrix(A, type = gpu_type)
  
  Q_prime <- qr.Q(qr(A))
  Wf <- function(v){
    v - (1 - sqrt(gamma)) * Q_prime %*% crossprod(Q_prime, v)
  }
  return(Wf)
}

Qf_temp <- function(v, Ue, Qf){
  Qfv <- Qf(v)
  Qfv - Ue %*% crossprod(Ue, Qfv)
}
