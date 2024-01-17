# dependencies:
library(parallel)
library(data.tree)
library(RcppEigen)
library(data.table)## visualize trees
library(igraph)
library(data.tree)
#' @importFrom Rdpack reprompt

n_cores <- detectCores()
# if there are less than 24 cores, it will be a local machine leave two cores for other tasks
if(n_cores > 1 && n_cores <= 24){
    n_cores <- n_cores - 1
}


#' Estimation of spectral transformation
#' 
#' Estimates the spectral transformation Q for spectral deconfounding by shrinking the leading singular values of the covariates.
#' @references
#'   \insertAllCited{}
#' @param X numerical covariates of class \code{matrix}
#' @param type type of deconfounding, one of 'trim', 'DDL_trim', 'pca', 'no_deconfounding'. 'trim' corresponds to the Trim transform \insertCite{Cevid2020SpectralModels}{SDForest}, 'DDL_trim' to the implementation of the Doubly debiased lasso \insertCite{Guo2022DoublyConfounding}{SDForest}, 'pca' to the PCA transformation\insertCite{Paul2008PreconditioningProblems}{SDForest} and 'no_deconfounding' to the Identity
#' @param trim_quantile quantile for Trim transform and DDL Trim transform, only needed for trim and DDL_trim
#' @param confounding_dim assumed confounding dimension, only needed for pca
#' @return Q of class \code{matrix}, the spectral transformation
#' @examples
#' X <- matrix(rnorm(100 * 20), nrow = 100)
#' get_Q(X, 'trim')
#' get_Q(X, 'DDL_trim')
#' get_Q(X, 'pca', q_hat = 5)
#' get_Q(X, 'no_deconfounding')
#' @export
get_Q <- function(X, type, trim_quantile = 0.5, confounding_dim = 0){
  # X: covariates
  # type: type of deconfounding
  modes <- c('trim' = 1, 'DDL_trim' = 2, 'pca' = 3, 'no_deconfounding' = 4)
  if(!(type %in% names(modes))) stop(paste("type must be one of:", paste(names(modes), collapse = ', ')))

  # number of observations
  n <- dim(X)[1]

  # calculate deconfounding matrix
  sv <- svd(X)
  tau <- quantile(sv$d, trim_quantile)
  D_tilde <- unlist(lapply(sv$d, FUN = function(x)min(x, tau))) / sv$d

  Q <- switch(modes[type], sv$u %*% diag(D_tilde) %*% t(sv$u), # trim
                          diag(n) - sv$u %*% diag(1 - D_tilde) %*% t(sv$u), # DDL_trim
                          { # pca
                              d_pca <- sv$d
                              if(confounding_dim <= 0) stop("the assumed confounding dimension must be larger than zero")
                              d_pca[1:confounding_dim] <- 0
                              sv$u %*% diag(d_pca) %*% t(sv$u)
                          },
                         diag(n)) # no_deconfounding
  return(Q)
}


##### SDTree functions #####

#' Spectral Deconfounded Tree
#' 
#' @param formula Object of class \code{formula} or describing the model to fit of the form \code{y ~ x1 + x2 + ...} where \code{y} is a numeric response and \code{x1, x2, ...} are vectors of covariates. Interactions are not supported.
#' @param data \code{data.frame} containing the variables in the model
#' @param x \code{matrix} of covariates
#' @param y \code{numeric} response
#' @param max_leaves maximum number of leaves for the grown tree
#' @param cp complexity parameter, minimum loss decrease to split a node. A split is only performed if the loss decrease is larger than \code{cp * initial_loss}, where \code{initial_loss} is the loss of the initial estimate using only a stump.
#' @param min_sample minimum number of observations per leaf
#' @param mtry number of randomly selected covariates to consider for a split, if \code{NULL} all covariates are available for each split
#' @param fast if \code{TRUE} only the optimal splitts in the new leaves are reevaluated, if \code{FALSE} all possible splitts are reevaluated
#' @param multicore if \code{TRUE} the optional splitts are evaluated in parallel (only available on unix systems)
#' @param mc.cores number of cores to use for parallel computing, if \code{NULL} all available cores - 1 are used
#' @param Q_type type of spectral transformation, see \code{\link{get_Q}}
#' @param trim_quantile quantile for Trim transform and DDL Trim transform, only needed for trim and DDL_trim, see \code{\link{get_Q}}
#' @param confounding_dim assumed confounding dimension, only needed for pca, see \code{\link{get_Q}}
#' @param Q spectral transformation, if \code{NULL} it is estimated using \code{\link{get_Q}}
#' @return Object of class \code{SDTree} containing
#' \item{predictions}{predictions for the training set}
#' \item{tree}{\code{data.tree} object containing the tree}
#' \item{var.names}{names of the covariates}
#' @examples
#' 
SDTree <- function(formula = NULL, data = NULL, x = NULL, y = NULL, max_leaves = 50, cp = 0.01, min_sample = 5, mtry = NULL, fast = TRUE,
                   multicore = F, mc.cores = NULL, Q_type = 'trim', trim_quantile = 0.5, confounding_dim = 0, Q = NULL){
  # X: covariates
  # Y: outcome
  # m: max number of splits
  # min_sample: min number of samples per leaf to split
  # cp: complexity parameter
  # mtry: if a number, only a random set of size mtry of the covariates are usable for a split
  input_data <- data.handler(formula = formula, data = data, x = x, y = y)
  X <- input_data$X
  Y <- input_data$Y

  # number of observations
  n <- dim(X)[1]
  # number of covariates
  p <- dim(X)[2]

  m <- max_leaves - 1
  # check validity of input
  if(n != length(Y)) stop('X and Y must have the same number of observations')
  if(m < 1) stop('max_leaves must be larger than 1')
  if(min_sample < 1) stop('min_sample must be larger than 0')
  if(cp < 0) stop('cp must be at least 0')
  if(!is.null(mtry) && mtry < 1) stop('mtry must be larger than 0')
  if(!is.null(mtry) && mtry > p) stop('mtry must be at most p')

  # estimate spectral transformation
  if(is.null(Q)){
    Q <- get_Q(X, Q_type, trim_quantile, confounding_dim)
  }else{
    if(!is.matrix(Q)) stop('Q must be a matrix')
    if(any(dim(Q) != n)) stop('Q must have dimension n x n')
  }

  # calculate first estimate
  E <- matrix(1, n, 1)
  E_tilde <- matrix(rowSums(Q))
  Y_tilde <- Q %*% Y

  # solve linear model
  c_hat <- fastLmPure(E_tilde, Y_tilde)$coefficients

  loss_start <- sum((Y_tilde - c_hat) ** 2) / n
  loss_temp <- loss_start

  # initialize tree
  tree <- Node$new(name = '1', value = c_hat, dloss = loss_start)

  # memory for optimal splits
  memory <- list(replicate(m + 1 , matrix(0, p, 2)))
  potential_splitts <- 1

  #available covariates
  if(!is.null(mtry)){
    len_p <- mtry
  }else {
    len_p <- p
  }

  for(i in 1:m){
    available_j <- sample(1:p, len_p)

    # iterate over all possible splits every time
    # for slow but slightly better solution
    if(!fast){
      potential_splitts <- 1:i
      to_small <- unlist(lapply(potential_splitts, function(x){sum(E[, x] == 1) < min_sample}))
      potential_splitts <- potential_splitts[!to_small]
    }

    #iterate over new to estimate splits
    for(branch in potential_splitts){
      best_splitts <- get_all_splitt(branch = branch, X = X, Y_tilde = Y_tilde, Q = Q, 
                                     n = n, n_branches = i+1, E = E, min_sample = min_sample, 
                                     multicore = multicore, mc.cores = mc.cores)
      best_splitts[, 1] <- loss_temp - best_splitts[, 1]
      memory[[branch]] <- best_splitts[, c(1, 3)]

      # if there is no available split fix loss decrease to zero
      if(dim(memory[[branch]])[1] == 0){
        memory[[branch]] <- matrix(0, p, 2)
      }
    }

    # find best split and its loss decrease
    Losses_dec <- unlist(lapply(memory, function(branch){branch[available_j, 1]}))
    loc <- which.max(Losses_dec)
    best_branch <- ceiling(loc / len_p)

    what_j <- loc %% len_p
    if(what_j == 0){
      what_j <- len_p
    }

    j <- available_j[what_j]
    s <- memory[[best_branch]][j, 2]

    # divide observations in leave
    index <- which(E[, best_branch] == 1)
    index_n_branches <- index[X[index, j] > s]

    # new indicator matrix
    E <- cbind(E, matrix(0, n, 1))
    E[index_n_branches, best_branch] <- 0
    E[index_n_branches, i+1] <- 1

    # calculate new level estimates
    E_tilde <- Q %*% E
    c_hat <- fastLmPure(E_tilde, Y_tilde)$coefficients

    # check if loss decrease is larger than minimum loss decrease
    # and if linear model could be estimated
    if(sum(is.na(c_hat)) > 0){
      warning('singulaer matrix QE, tree might be to large, consider increasing cp')
      break
    }
    if(Losses_dec[loc] <= cp * loss_start){
      break
    }

    # select leave to split
    if(tree$height == 1){
      leave <- tree
    }else{
      leaves <- tree$leaves
      leave <- leaves[[which(tree$Get('name', filterFun = isLeaf) == best_branch)]]
    }

    # save split rule
    leave$j <- j
    leave$s <- s

    # add new leaves
    leave$AddChild(best_branch, value = 0, dloss = Losses_dec[loc], decision = 'no')
    leave$AddChild(i + 1, value = 0, dloss = Losses_dec[loc], decision = 'yes')

    # add estimates to tree leaves
    for(l in tree$leaves){
      l$value <- c_hat[as.numeric(l$name)]
    }

    # new temporary loss
    loss_temp <- loss(Y_tilde, E_tilde %*% c_hat)

    # the two new partitions need to be checked for optimal splits in next iteration
    potential_splitts <- c(best_branch, i + 1)

    # a partition with less than min_sample observations or unique samples are not available for further splits
    to_small <- unlist(lapply(potential_splitts, function(x){(sum(E[, x] == 1) < min_sample)| length(unique(X[which(E[, x] == 1), 1])) == 1}))
    if(sum(to_small) > 0){
      for(el in potential_splitts[to_small]){
        # to small partitions cannot decrease the loss
        memory[[el]] <- matrix(0, p, 2)
      }
      potential_splitts <- potential_splitts[!to_small]
    }
  }

  # print warning if maximum splitts was reached, one might want to increase m
  if(i == m){
    warning('maximum number of iterations was reached, consider increasing m!')
  }

  # predict the test set
  f_X_hat <- predict_outsample(tree, X)

  # labels for the nodes
  tree$Do(splitt_names, filterFun = isNotLeaf, var_names = colnames(X))
  tree$Do(leave_names, filterFun = isLeaf)

  res <- list(predictions = f_X_hat, tree = tree, var.names = colnames(X))
  class(res) <- 'SDTree'
  return(res)
}

cv.SDTree <- function(formula = NULL, data = NULL, x = NULL, y = NULL, max_leaves = 50, cp = 0, 
                      min_sample = 5, fast = TRUE, multicore = F, Q_type = 'trim', trim_quantile = 0.5, 
                      confounding_dim = 0, mc.cores = NULL, n_cv = 3, cp.seq = NULL){

  input_data <- data.handler(formula = formula, data = data, x = x, y = y)
  X <- input_data$X
  Y <- input_data$Y

  p <- ncol(X)
  n <- nrow(X)

  m <- max_leaves - 1
  # check validity of input
  if(n != length(Y)) stop('X and Y must have the same number of observations')
  if(m < 1) stop('max_leaves must be larger than 1')
  if(min_sample < 1) stop('min_sample must be larger than 0')
  if(cp < 0) stop('cp must be at least 0')
  if(!is.null(cp.seq) && (any(cp.seq < 0) | any(cp.seq > 1))) stop('cp.seq must be between 0 and 1')
  if(n_cv < 2 | n_cv > n) stop('n_cv must be at least 2 and smaller than n')
  if(n_cv > 5) warning('if to many folds are used, Q_validation might differ to much from Q_trainig, consider using less folds')

  # estimate spectral transformation
  Q <- get_Q(X, Q_type, trim_quantile, confounding_dim)
  
  # estimating initial loss with only a stump
  # to map optimal minimal loss decrease to a cp value
  E <- matrix(1, n, 1)
  E_tilde <- matrix(rowSums(Q))
  Y_tilde <- Q %*% Y
  # solve linear model
  c_hat <- fastLmPure(E_tilde, Y_tilde)$coefficients
  loss_start <- sum((Y_tilde - c_hat) ** 2) / n

  # validation set size
  len_test <- floor(n / n_cv)
  # validation sets
  test_ind <- lapply(1:n_cv, function(x)1:len_test + (x - 1) * len_test)

  # sequence of minimum loss decrease to compare with cv
  if(is.null(cp.seq)){
    cp.seq <- seq(0, 6, 0.002)
  }
  t_seq <- cp.seq * loss_start

  # estimate performance for every validation set
  perf <- lapply(test_ind, function(cv_ind){
    # calculate Trim transform
    Q_cv <- get_Q(X[cv_ind, ], Q_type, trim_quantile, confounding_dim)

    # deconfound X and Y
    X_train <- X[-cv_ind, ]
    Y_train <- Y[-cv_ind]

    X_cv <- X[cv_ind, ]
    Y_cv <- Y[cv_ind]
    
    # estimate tree with the training set
    suppressWarnings({
    res <- SDTree(x = X_train, y = Y_train, max_leaves = max_leaves, cp = cp, min_sample = min_sample,
                  Q_type = Q_type, trim_quantile = trim_quantile, confounding_dim = confounding_dim,
                  multicore = multicore, fast = fast, mc.cores = mc.cores)
    })

    # validation performance if we prune with the different ts
    if(multicore){
      if(!is.null(mc.cores)){
        n_cores <- mc.cores
      }
      perf <- mclapply(t_seq, function(t) pruned_loss(res$tree, X_cv, Y_cv, Q_cv, t), mc.cores = n_cores)
    }else{
      perf <- lapply(t_seq, function(t) pruned_loss(res$tree, X_cv, Y_cv, Q_cv, t))
    }
    
    return(perf)
  })
  
  # collect performance for different min loss decreases
  perf <- matrix(unlist(perf), ncol = n_cv, byrow = FALSE)
  
  cp.table <- matrix(c(t_seq / loss_start, apply(perf, 1, mean), apply(perf, 1, sd)), ncol = 3, byrow = FALSE)
  colnames(cp.table) <- c('cp', 'SDLoss mean', 'SDLoss sd')

  loss_unique <- unique(cp.table[, 2])
  cp.table <- lapply(loss_unique, function(loss){
      idx <- which(cp.table[, 2] == loss)
      cp.table[idx[length(idx)], ]
  })
  cp.table <- do.call(rbind, cp.table)

  cp.min <- cp.table[which.min(cp.table[, 2]), 1]

  res <- list(cp.min = cp.min, cp.table = cp.table)
  return(res)
}

predict.SDTree <- function(object, newdata){
  # predict function for the spectral deconfounded tree
  if(!is.data.frame(newdata)) stop('newdata must be a data.frame')
  if(!all(object$var.names %in% names(newdata))) stop('newdata must contain all covariates used for training')

  X <- newdata[, object$var.names]
  return(predict_outsample(object$tree, X))
}

print.SDTree <- function(x){
  # print function for the spectral deconfounded tree
  print(x$tree, 'value', 's', 'j', 'label', 'decision')
}

plot.SDTree <- function(x){
  # plot function for the spectral deconfounded tree
  SetEdgeStyle(x$tree, label = function(x) {x$decision})
  SetNodeStyle(x$tree, label = function(x) {x$label})
  plot(x$tree)
}

#### SDForest functions ####

# TODO: add variable importance
# TODO: add oob error
SDForest <- function(formula = NULL, data = NULL, x = NULL, y = NULL, nTree = 100, max_leaves = 50, 
                     cp = 0, min_sample = 5, mtry = NULL, multicore = F, mc.cores = NULL, 
                     Q_type = 'trim', trim_quantile = 0.5, confounding_dim = 0, Q = NULL){

  input_data <- data.handler(formula = formula, data = data, x = x, y = y)
  X <- input_data$X
  Y <- input_data$Y

  # number of observations
  n <- nrow(X)

  if(n != length(Y)) stop('X and Y must have the same number of observations')
  if(!is.null(mtry) && mtry < 1) stop('mtry must be larger than 0')
  if(!is.null(mtry) && mtry > p) stop('mtry must be at most p')

  # estimate spectral transformation
  if(is.null(Q)){
    Q <- get_Q(X, Q_type, trim_quantile, confounding_dim)
  }else{
    if(!is.matrix(Q)) stop('Q must be a matrix')
    if(any(dim(Q) != n)) stop('Q must have dimension n x n')
  }

  # mtry
  if(is.null(mtry)){
    mtry <- floor(sqrt(ncol(X)))
  }
  
  # bootstrap samples
  ind <- lapply(1:nTree, function(x)sample(1:n, n, replace = T))
  
  suppressWarnings({
  # estimating all the trees
  if(multicore){
    if(!is.null(mc.cores)){
      n_cores <- mc.cores
    }
    res <- mclapply(ind, function(i)SDTree(x = X[i, ], y = Y[i], max_leaves = max_leaves, cp = cp, 
                                           min_sample = min_sample, Q = Q, mtry = mtry),
                    mc.cores = n_cores)
  }else{
    res <- lapply(ind, function(i)SDTree(x = X[i, ], y = Y[i], max_leaves = max_leaves, cp = cp, 
                                           min_sample = min_sample, Q = Q, mtry = mtry))
  }
  })
  
  # predict with all trees
  pred <- do.call(cbind, lapply(res, function(x){predict_outsample(x$tree, X)}))
  
  # use mean over trees as final prediction
  f_X_hat <- rowMeans(pred)
  res <- list(predictions = f_X_hat, forest = res, var.names = colnames(X))
  class(res) <- 'SDforest'
  return(res)
}

predict.SDforest <- function(object, newdata){
  # predict function for the spectral deconfounded random forest
  # using the mean over all trees as the prediction
  # check data type
  if(!is.data.frame(newdata)) stop('newdata must be a data.frame')
  if(!all(object$var.names %in% names(newdata))) stop('newdata must contain all covariates used for training')

  X <- newdata[, object$var.names]

  pred <- do.call(cbind, lapply(object[[2]], function(x){predict_outsample(x$tree, X)}))
  return(rowMeans(pred))
}

#TODO: print function

#### Utility functions ####

#helper functions to lable nodes
splitt_names <- function(node, var_names = NULL){
  if(is.null(var_names)){
    node$label <- paste('X', node$j, ' <= ', round(node$s, 2), sep = '')
  }else{
    node$label <- paste(var_names[node$j], ' <= ', round(node$s, 2), sep = '')
  }
}

leave_names <- function(node){
    new_name <- as.character(round(node$value, 1))
    if(new_name %in% node$Get('name', filterFun = isLeaf)){
        new_name <- paste(new_name, '')
    }
    node$label <- new_name
}

evaluate_splitt <- function(branch, j, s, index, X, Y_tilde, Q, n, n_branches, E, min_sample){
  # evaluate a split at partition branch on covariate j at the splitpoint s
  # index: index of observations in branch
  
  # dividing observation in branch
  index_n_branches <- index[X[index, j] > s]
  index_branch <- index[X[index, j] <= s]
  
  # check wether this split resolves in two reasnable partitions
  if(length(index_branch) < min_sample | length(index_n_branches) < min_sample){
    # remove no longer needed objects from memory
    rm(Q, index, index_n_branches)
    return(list('loss' = Inf, j = j, s = s))
  }
  
  # new indicator matrix
  E <- cbind(E, matrix(0, n, 1))
  E[index_n_branches, branch] <- 0
  E[index_n_branches, n_branches] <- 1
  
  # new level estimates
  E_tilde <- Q %*% E
  c_hat <- fastLmPure(E_tilde, Y_tilde)$coefficients

  if(any(is.na(c_hat))){
    # linear model could not be estimated
    loss <- Inf
  }else{
    # resulting new spectral loss
    loss <- loss(Y_tilde, E_tilde %*% c_hat)
  }

  # remove no longer needed objects from memory
  rm(Q, Y_tilde, E, E_tilde, c_hat, index, index_n_branches)
  return(list('loss' = loss, j = j, s = s))
}

get_all_splitt <- function(branch, X, Y_tilde, Q, n, n_branches, E, min_sample, multicore, mc.cores = NULL){
  # finds the best splitts for every covariate in branch
  # returns the best splitpoint for every covariate and the resulting loss decrease

  if(!is.null(mc.cores)){
    n_cores <- mc.cores
  }

  # observations belonging to branch
  index <- which(E[, branch] == 1)

  # all possible split points
  s <- find_s(X[index, ])
  # itetator for the splits
  iter <- 1:length(s) - 1

  # evaluate all the relevant splitts
  if(multicore){
      res <- mclapply(iter, function(x)evaluate_splitt(branch = branch, j = floor(x / dim(s)[1] + 1), 
              s = s[x + 1], index = index, X = X, Y_tilde = Y_tilde, Q = Q, n = n, n_branches = n_branches, 
              E = E, min_sample = min_sample), mc.cores = n_cores)
  }else{
    res <- lapply(iter, function(x)evaluate_splitt(branch = branch, j = floor(x / dim(s)[1] + 1), 
              s = s[x + 1], index = index, X = X, Y_tilde = Y_tilde, Q = Q, n = n, n_branches = n_branches, 
              E = E, min_sample = min_sample))
  }

  # choose minimum split for every covariate
  res <- rbindlist(res)
  res_min <- lapply(1:ncol(X), function(x){res_temp <- res[j == x, ]
                                           res_temp[which.min(loss), ]})
  res_min <- matrix(unlist(res_min), ncol = 3, byrow = T)
  return(res_min)
}

find_s <- function(X){
  # finds all the reasnable splitting points in a data matrix

  if(is.null(dim(X))){
    X <- matrix(X, ncol = 1)
  }

  X_sort <- apply(X, 2, sort)

  # find middlepoints between observed x values
  s <- X_sort[-dim(X)[1], ] + diff(X_sort)/2

  # for runtime reasons use only every 4th middlepoint if we have more than 400 observations
  # and only every second if we have more than 200 observations
  if(dim(s)[1] > 400){
    s <- s[seq(1, dim(s)[1], 4), ]
  }else if (dim(s)[1] > 200) {
    s <- s[seq(1, dim(s)[1], 2), ]
  }
  return(s)
}

traverse_tree <- function(tree, x){
  # traverse the tree using the splitting rules and 
  # returns point estimate for f(x)
  if(tree$isLeaf){
    return(tree$value)
  }
  if(x[tree$j] <= tree$s){
    return(traverse_tree(tree$children[[1]], x))
  }else {
    return(traverse_tree(tree$children[[2]], x))
  }
}

predict_outsample <- function(tree, X){
  # predict for every observation in X f(x)
  # using the splitting rules from the tree
  return(apply(X, 1, function(x)traverse_tree(tree, x)))
}

loss <- function(Y, f_X){
  # MSE
  return(sum((Y - f_X)^2) / length(Y))
}

pruned_loss <- function(tree, X_val, Y_val, Q_val, t){
  # funtion to prune tree using the minimum loss decrease t
  # and return spectral loss on the validation set
  
  tree_t <- Clone(tree)
  
  # prune tree
  Prune(tree_t, function(x) x$dloss > t)
  
  # predict on test set
  f_X_hat_val <- predict_outsample(tree_t, X_val)
  
  # return spectral loss
  return(sum((Q_val %*% Y_val - Q_val %*% f_X_hat_val) ** 2) / length(Y_val))
}

data.handler <- function(formula = NULL, data = NULL, x = NULL, y = NULL){
  if(is.null(formula)){
    if(is.null(x) | is.null(y)){
      stop("Error: Either data or x and y is required.")
    }else {
      x <- apply(x, 2, function(x){if (is.character(x)) as.numeric(factor(x))
                                    else if(!is.numeric(x))  as.numeric(x)
                                    else x})
      if (!is.numeric(y)) stop("Error: y must be numeric. Only regression is supported at the moment.")
      return(list(X = as.matrix(x), Y = y))
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
      return(list(X = X, Y = Y))
    }
  }
}