# dependencies:
#library(parallel)
#library(doParallel)
#library(foreach)
#library(data.tree)
#library(RcppEigen)
#library(data.table)## visualize trees
#library(igraph)
#library(data.tree)
#' @importFrom Rdpack reprompt

n_cores <- parallel::detectCores()
 # if there are less than 24 cores, it will be a local machine leave two cores for other tasks
if(n_cores > 1 && n_cores <= 24){
    n_cores <- n_cores - 1
}


#' Estimation of spectral transformation
#' 
#' Estimates the spectral transformation Q for spectral deconfounding by shrinking the leading singular values of the covariates.
#' @references
#'   \insertAllCited{}
#' @author Markus Ulmer
#' @param X Numerical covariates of class \code{matrix}.
#' @param type Type of deconfounding, one of 'trim', 'DDL_trim', 'pca', 'no_deconfounding'. 'trim' corresponds to the Trim transform \insertCite{Cevid2020SpectralModels}{SDForest}, 'DDL_trim' to the implementation of the Doubly debiased lasso \insertCite{Guo2022DoublyConfounding}{SDForest}, 'pca' to the PCA transformation\insertCite{Paul2008PreconditioningProblems}{SDForest} and 'no_deconfounding' to the Identity.
#' @param trim_quantile Quantile for Trim transform and DDL Trim transform, only needed for trim and DDL_trim.
#' @param confounding_dim Assumed confounding dimension, only needed for pca.
#' @return Q of class \code{matrix}, the spectral transformation matrix.
#' @examples
#' X <- matrix(rnorm(100 * 20), nrow = 100)
#' get_Q(X, 'trim')
#' get_Q(X, 'DDL_trim')
#' get_Q(X, 'pca', q_hat = 5)
#' get_Q(X, 'no_deconfounding')
#' @export
get_Q <- function(X, type, trim_quantile = 0.5, confounding_dim = 0){
  svd_error <- function(X, f = 1, count = 1){
    tryCatch({
      svd(X * f)
    }, error = function(e) {
        warning(paste(e, ':X multipied by number close to 1'))
        if(count > 5) stop('svd did not converge')
        return(svd_error(X, 1 + 0.0000000000000001 * 10 ^ count, count + 1))})
  }

  # X: covariates
  # type: type of deconfounding
  modes <- c('trim' = 1, 'pca' = 2, 'no_deconfounding' = 3)
  if(!(type %in% names(modes))) stop(paste("type must be one of:", paste(names(modes), collapse = ', ')))

  # number of observations
  n <- dim(X)[1]

  # calculate deconfounding matrix
  sv <- svd_error(X)

  tau <- quantile(sv$d, trim_quantile)
  D_tilde <- unlist(lapply(sv$d, FUN = function(x)min(x, tau))) / sv$d
  D_tilde[is.na(D_tilde)] <- 1

  Q <- switch(modes[type], diag(n) - sv$u %*% diag(1 - D_tilde) %*% t(sv$u), # DDL_trim
                          { # pca
                              d_pca <- rep(1, length(sv$d))
                              if(confounding_dim <= 0) stop("the assumed confounding dimension must be larger than zero")
                              d_pca[1:confounding_dim] <- 0
                              diag(n) - sv$u %*% diag(1 - d_pca) %*% t(sv$u)
                          },
                         diag(n)
                         ) # no_deconfounding
  return(Q)
}

condDependence <- function(object, j, X = NULL, multicore = F, mc.cores = NULL){
  if(is.character(j)){
    j <- which(names(X) == j)
  }

  if(is.null(X)){
    X <- data.frame(object$X)
    if(is.null(X)) stop('X must be provided if it is not part of the object')
  }

  if(!is.numeric(j)) stop('j must be a numeric or character')
  if(j > ncol(X)) stop('j must be smaller than p')
  if(j < 1) stop('j must be larger than 0')
  if(any(is.na(X))) stop('X must not contain missing values')

  #x_seq <- seq(quantile(X[, j], 0.05), quantile(X[, j], 0.95), length.out = 100)
  x_seq <- seq(min(X[, j]), max(X[, j]), length.out = 100)

  if(multicore){
    if(!is.null(mc.cores)){
      n_cores <- mc.cores
    }
    preds <- parallel::mclapply(x_seq, function(x){
      X_new <- X
      X_new[, j] <- x
      pred <- predict(object, newdata = X_new)
      return(pred)
    }, mc.cores = n_cores)
  }else{
    preds <- pbapply::pblapply(x_seq, function(x){
      X_new <- X
      X_new[, j] <- x
      pred <- predict(object, newdata = X_new)
      return(pred)
    })
  }
  preds <- do.call(rbind, preds)
  preds_mean <- rowMeans(preds)

  res <- list(preds_mean = preds_mean, x_seq = x_seq, preds = preds, j = j, xj = X[, j])
  class(res) <- 'condDependence'
  return(res)
}

plot.condDependence <- function(object, n_examples = 19){
  ggdep <- ggplot2::ggplot() + ggplot2::theme_bw()
  preds <- object$preds
  x_seq <- object$x_seq
  
  sample_examples <- sample(1:dim(preds)[2], n_examples)
  for(i in sample_examples){
      pred_data <- data.frame(x = x_seq, y = preds[, i])
      ggdep <- ggdep + ggplot2::geom_line(data = pred_data, ggplot2::aes(x = x, y = y), col = 'grey')
  }

  ggdep <- ggdep + ggplot2::geom_line(data = data.frame(x = x_seq, y = object$preds_mean), 
                    ggplot2::aes(x = x, y = y), col = '#08cbba', linewidth = 1.5)
  ggdep <- ggdep + ggplot2::geom_point(data = data.frame(x = object$xj, y = min(preds[, sample_examples])), 
                    ggplot2::aes(x = x, y = y), col = 'black', size = 1,shape = 108)
  ggdep <- ggdep + ggplot2::ylab('f(x)') + ggplot2::ggtitle('Conditional dependence')
  ggdep <- ggdep + ggplot2::xlim(quantile(object$xj, 0.05), quantile(object$xj, 0.95))
  if(is.character(object$j)){
    ggdep <- ggdep + ggplot2::xlab(object$j)
  }else{
    ggdep <- ggdep + ggplot2::xlab(paste('x', object$j, sep = ''))
  }
  ggdep
}



##### SDTree functions #####

#' Spectral Deconfounded Tree
#' 
#' Estimates a tree using spectral deconfounding. # TODO: add more details
#' @references
#'  \insertAllCited{}
#' @author Markus Ulmer
#' @param formula Object of class \code{formula} or describing the model to fit of the form \code{y ~ x1 + x2 + ...} where \code{y} is a numeric response and \code{x1, x2, ...} are vectors of covariates. Interactions are not supported.
#' @param data Training data of class \code{data.frame} containing the variables in the model.
#' @param x Predictor data, alternative to \code{formula} and \code{data}.
#' @param y Response vector, alternative to \code{formula} and \code{data}.
#' @param max_leaves Maximum number of leaves for the grown tree.
#' @param cp Complexity parameter, minimum loss decrease to split a node. A split is only performed if the loss decrease is larger than \code{cp * initial_loss}, where \code{initial_loss} is the loss of the initial estimate using only a stump.
#' @param min_sample Minimum number of observations per leaf. A split is only performed if both resulting leaves have at least \code{min_sample} observations.
#' @param mtry Number of randomly selected covariates to consider for a split, if \code{NULL} all covariates are available for each split.
#' @param fast If \code{TRUE}, only the optimal splitts in the new leaves are evaluated and the previously optimal splitts and their potential loss-decrease are reused. If \code{FALSE} all possible splitts in all the leaves are reevaluated after every split.
#' @param multicore If \code{TRUE} the optional splitts are evaluated in parallel (only available on unix systems).
#' @param mc.cores Number of cores to use for parallel computing, if \code{NULL} all available cores - 1 are used.
#' @param Q_type Type of deconfounding, one of 'trim', 'DDL_trim', 'pca', 'no_deconfounding'. 'trim' corresponds to the Trim transform \insertCite{Cevid2020SpectralModels}{SDForest}, 'DDL_trim' to the implementation of the Doubly debiased lasso \insertCite{Guo2022DoublyConfounding}{SDForest}, 'pca' to the PCA transformation\insertCite{Paul2008PreconditioningProblems}{SDForest} and 'no_deconfounding' to the Identity. See \code{\link{get_Q}}.
#' @param trim_quantile Quantile for Trim transform and DDL Trim transform, only needed for trim and DDL_trim, see \code{\link{get_Q}}.
#' @param confounding_dim Assumed confounding dimension, only needed for pca, see \code{\link{get_Q}}.
#' @param Q Spectral transformation, if \code{NULL} it is internally estimated using \code{\link{get_Q}}.
#' @return Object of class \code{SDTree} containing
#' \item{predictions}{Predictions for the training set.}
#' \item{tree}{The estimated tree of type \code{data.tree} \insertCite{data.tree}{SDForest}. The tree contains the information about all the splits and the resulting estimates.}
#' \item{var_names}{Names of the covariates in the training data.}
#' @examples
#' # TODO: add example
#' @export
SDTree <- function(formula = NULL, data = NULL, x = NULL, y = NULL, max_leaves = 50, cp = 0.01, min_sample = 5, mtry = NULL, fast = TRUE,
                   Q_type = 'trim', trim_quantile = 0.5, confounding_dim = 0, Q = NULL){
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
  if(n < 2 * min_sample) stop('n must be at least 2 * min_sample')

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

  u_start <- E_tilde / sqrt(sum(E_tilde ** 2))
  Q_temp <- Q - u_start %*% (t(u_start) %*% Q)

  Y_tilde <- Q %*% Y
  
  # solve linear model
  c_hat <- RcppEigen::fastLmPure(E_tilde, Y_tilde)$coefficients

  #c_hat <- lm.fit(E_tilde, Y_tilde)$coefficients

  loss_start <- sum((Y_tilde - c_hat) ** 2) / n
  loss_temp <- loss_start

  # initialize tree
  tree <- data.tree::Node$new(name = '1', value = c_hat, dloss = loss_start, cp = 10)

  # memory for optimal splits
  memory <- list(replicate(m + 1 , matrix(0, p, 4)))
  potential_splitts <- 1

  # variable importance
  var_imp <- rep(0, p)
  names(var_imp) <- colnames(X)

  after_mtry <- 0

  for(i in 1:m){
    # iterate over all possible splits every time
    # for slow but slightly better solution
    if(!fast){
      potential_splitts <- 1:i
      to_small <- unlist(lapply(potential_splitts, function(x){sum(E[, x] == 1) < min_sample*2}))
      potential_splitts <- potential_splitts[!to_small]
    }
    #iterate over new to estimate splits
    for(branch in potential_splitts){
      best_splitts <- get_all_splitt(branch = branch, X = X, Y_tilde = Y_tilde, Q_temp = Q_temp, 
                                     n = n, min_sample = min_sample, p = p, E = E)
      memory[[branch]] <- best_splitts
    }

    if(i > after_mtry && !is.null(mtry)){
      Losses_dec <- lapply(memory, function(branch){
        branch[sample(1:p, mtry), ]})
      Losses_dec <- do.call(rbind, Losses_dec)
    }else {
       Losses_dec <- do.call(rbind, memory)
    }

    if(max(Losses_dec[, 1]) <= 0){
      break
    }
    loc <- which.max(Losses_dec[, 1])
    best_branch <- Losses_dec[loc, 4]
    j <- Losses_dec[loc, 2]
    s <- Losses_dec[loc, 3]

    # divide observations in leave
    index <- which(E[, best_branch] == 1)
    index_n_branches <- index[X[index, j] > s]


    # new indicator matrix
    E <- cbind(E, matrix(0, n, 1))
    E[index_n_branches, best_branch] <- 0
    E[index_n_branches, i+1] <- 1

    E_tilde_branch <- E_tilde[, best_branch]
    E_tilde[, best_branch] <- Q %*% E[, best_branch]
    E_tilde <- cbind(E_tilde, E_tilde_branch - E_tilde[, best_branch])

    c_hat <- RcppEigen::fastLmPure(E_tilde, Y_tilde)$coefficients

    u_next_prime <- Q_temp %*% E[, i + 1]
    u_next <- u_next_prime / sqrt(sum(u_next_prime ** 2))

    Q_temp <- Q_temp - u_next %*% (t(u_next) %*% Q)

    # check if loss decrease is larger than minimum loss decrease
    # and if linear model could be estimated
    if(sum(is.na(c_hat)) > 0){
      warning('singulaer matrix QE, tree might be to large, consider increasing cp')
      break
    }

    loss_dec <- loss_temp - loss(Y_tilde, E_tilde %*% c_hat)
    loss_temp <- loss_temp - loss_dec

    if(loss_dec <= cp * loss_start){
      break
    }
    # add loss decrease to variable importance
    var_imp[j] <- var_imp[j] + loss_dec

    # select leave to split
    if(tree$height == 1){
      leave <- tree
    }else{
      leaves <- tree$leaves
      leave <- leaves[[which(tree$Get('name', filterFun = data.tree::isLeaf) == best_branch)]]
    }

    # save split rule
    leave$j <- j
    leave$s <- s

    leave$res_dloss <- loss_dec

    # add new leaves
    leave$AddChild(best_branch, value = 0, dloss = loss_dec, cp = loss_dec / loss_start, decision = 'no', n_samples = sum(E[, best_branch] == 1))
    leave$AddChild(i + 1, value = 0, dloss = loss_dec, cp = loss_dec / loss_start, decision = 'yes', n_samples = sum(E[, i + 1] == 1))

    # add estimates to tree leaves
    for(l in tree$leaves){
      l$value <- c_hat[as.numeric(l$name)]
    }

    # the two new partitions need to be checked for optimal splits in next iteration
    potential_splitts <- c(best_branch, i + 1)

    # a partition with less than min_sample observations or unique samples are not available for further splits
    to_small <- unlist(lapply(potential_splitts, function(x){length(unique(X[which(E[, x] == 1), 1])) < min_sample * 2}))
    if(sum(to_small) > 0){
      for(el in potential_splitts[to_small]){
        # to small partitions cannot decrease the loss
        memory[[el]] <- matrix(0, p, 4)
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

  var_names <- colnames(data.frame(X))
  names(var_imp) <- var_names

  # labels for the nodes
  tree$Do(splitt_names, filterFun = data.tree::isNotLeaf, var_names = var_names)
  tree$Do(leave_names, filterFun = data.tree::isLeaf)

  res <- list(predictions = f_X_hat, tree = tree, var_names = var_names, var_importance = var_imp)
  class(res) <- 'SDTree'
  return(res)
}

#' Cross-validation for the Spectral Deconfounded Tree
#' 
#' Estimates the optimal complexity parameter for the spectral deconfounded tree using cross-validation. Q is estimated for each training set and validation set separately to ensure independence of the validation set.
#' @references
#' \insertAllCited{}
#' @author Markus Ulmer
#' @param formula Object of class \code{formula} or describing the model to fit of the form \code{y ~ x1 + x2 + ...} where \code{y} is a numeric response and \code{x1, x2, ...} are vectors of covariates. Interactions are not supported.
#' @param data Training data of class \code{data.frame} containing the variables in the model.
#' @param x Predictor data, alternative to \code{formula} and \code{data}.
#' @param y Response vector, alternative to \code{formula} and \code{data}.
#' @param max_leaves Maximum number of leaves for the grown tree.
#' @param min_sample Minimum number of observations per leaf. A split is only performed if both resulting leaves have at least \code{min_sample} observations.
#' @param mtry Number of randomly selected covariates to consider for a split, if \code{NULL} all covariates are available for each split.
#' @param fast If \code{TRUE}, only the optimal splitts in the new leaves are evaluated and the previously optimal splitts and their potential loss-decrease are reused. If \code{FALSE} all possible splitts in all the leaves are reevaluated after every split.
#' @param multicore If \code{TRUE} the optional splitts are evaluated in parallel (only available on unix systems).
#' @param mc.cores Number of cores to use for parallel computing, if \code{NULL} all available cores - 1 are used.
#' @param Q_type Type of deconfounding, one of 'trim', 'DDL_trim', 'pca', 'no_deconfounding'. 'trim' corresponds to the Trim transform \insertCite{Cevid2020SpectralModels}{SDForest}, 'DDL_trim' to the implementation of the Doubly debiased lasso \insertCite{Guo2022DoublyConfounding}{SDForest}, 'pca' to the PCA transformation\insertCite{Paul2008PreconditioningProblems}{SDForest} and 'no_deconfounding' to the Identity. See \code{\link{get_Q}}.
#' @param trim_quantile Quantile for Trim transform and DDL Trim transform, only needed for trim and DDL_trim, see \code{\link{get_Q}}.
#' @param confounding_dim Assumed confounding dimension, only needed for pca, see \code{\link{get_Q}}.
#' @param n_cv Number of folds for cross-validation. It is recommended to not use more than 5 folds if the number of covariates is larger than the number of observations. In this case the spectral transformation could differ to much if the validation data is substantially smaller than the training data.
#' @param cp_seq Sequence of complexity parameters cp to compare using cross-validation, if \code{NULL} a sequence from 0 to 0.6 with stepsize 0.002 is used.
#' @return A list containing
#' \item{cp_min}{The optimal complexity parameter.}
#' \item{cp_table}{A table containing the complexity parameter, the mean and the standard deviation of the loss on the validation sets for the complexity parameters. If multiple complexity parameters result in the same loss, only the one with the largest complexity parameter is shown.}
#' @examples
#' # TODO: add example
#' @seealso \code{\link{SDTree}}
#' @export
cv.SDTree <- function(formula = NULL, data = NULL, x = NULL, y = NULL, max_leaves = 50, 
                      min_sample = 5, fast = TRUE, multicore = F, Q_type = 'trim', trim_quantile = 0.5, 
                      confounding_dim = 0, mc.cores = NULL, n_cv = 3, cp_seq = NULL){

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
  if(!is.null(cp_seq) && (any(cp_seq < 0) | any(cp_seq > 1))) stop('cp.seq must be between 0 and 1')
  if(n_cv < 2 | n_cv > n) stop('n_cv must be at least 2 and smaller than n')
  if(n_cv > 5 & n < p) warning('if n < p and to many folds are used, Q_validation might differ to much from Q_trainig, consider using less folds')

  # estimate spectral transformation
  Q <- get_Q(X, Q_type, trim_quantile, confounding_dim)
  
  # estimating initial loss with only a stump
  # to map optimal minimal loss decrease to a cp value
  E <- matrix(1, n, 1)
  E_tilde <- matrix(rowSums(Q))
  Y_tilde <- Q %*% Y

  # solve linear model
  c_hat <- RcppEigen::fastLmPure(E_tilde, Y_tilde)$coefficients
  loss_start <- sum((Y_tilde - c_hat) ** 2) / n

  # validation set size
  len_test <- floor(n / n_cv)
  # validation sets
  test_ind <- lapply(1:n_cv, function(x)1:len_test + (x - 1) * len_test)

  # sequence of minimum loss decrease to compare with cv
  if(is.null(cp_seq)){
    cp_seq <- seq(0, 0.6, 0.002)
  }
  t_seq <- cp_seq * loss_start

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
    res <- SDTree(x = X_train, y = Y_train, max_leaves = max_leaves, cp = 0, min_sample = min_sample,
                  Q_type = Q_type, trim_quantile = trim_quantile, confounding_dim = confounding_dim)
    })

    # validation performance if we prune with the different ts
    if(multicore){
      if(!is.null(mc.cores)){
        n_cores <- mc.cores
      }
      perf <- mclapply(t_seq, function(t) pruned_loss(res$tree, X_cv, Y_cv, Q_cv, t), mc.cores = n_cores, mc.preschedule = FALSE)
    }else{
      perf <- lapply(t_seq, function(t) pruned_loss(res$tree, X_cv, Y_cv, Q_cv, t))
    }
    
    return(perf)
  })
  
  # collect performance for different min loss decreases
  perf <- matrix(unlist(perf), ncol = n_cv, byrow = FALSE)
  
  cp_table <- matrix(c(t_seq / loss_start, apply(perf, 1, mean), apply(perf, 1, sd)), ncol = 3, byrow = FALSE)
  colnames(cp_table) <- c('cp', 'SDLoss mean', 'SDLoss sd')

  loss_unique <- unique(cp_table[, 2])
  cp_table <- lapply(loss_unique, function(loss){
      idx <- which(cp_table[, 2] == loss)
      cp_table[idx[length(idx)], ]
  })
  cp_table <- do.call(rbind, cp_table)

  cp_min <- cp_table[which.min(cp_table[, 2]), 1]

  res <- list(cp_min = cp_min, cp_table = cp_table)
  return(res)
}

#' Predictions for the Spectral Deconfounded Tree
#' 
#' Predicts the response for new data using a fitted spectral deconfounded tree.
#' @references
#' \insertAllCited{}
#' @author Markus Ulmer
#' @param object Fitted object of class \code{SDTree}.
#' @param newdata New test data of class \code{data.frame} containing the covariates for which to predict the response.
#' @return A vector of predictions for the new data.
#' @examples
#' # TODO: add example
#' @seealso \code{\link{SDTree}}
#' @export
predict.SDTree <- function(object, newdata){
  # predict function for the spectral deconfounded tree
  if(!is.data.frame(newdata)) stop('newdata must be a data.frame')
  if(!all(object$var_names %in% names(newdata))) stop('newdata must contain all covariates used for training')

  X <- newdata[, object$var_names]
  if(any(is.na(X))) stop('X must not contain missing values')
  return(predict_outsample(object$tree, X))
}

#' Print Spectral Deconfounded Tree
#' 
#' Print contents of the spectral deconfounded tree.
#' @references
#' \insertAllCited{}
#' @author Markus Ulmer
#' @param object Fitted object of class \code{SDTree}.
#' @seealso \code{\link{SDTree}}
print.SDTree <- function(object){
  # print function for the spectral deconfounded tree
  print(object$tree, 'value', 's', 'j', 'label', 'decision', 'n_samples')
}

#' Plot Spectral Deconfounded Tree
#' 
#' Plot the spectral deconfounded tree.
#' @references
#' \insertAllCited{}
#' @author Markus Ulmer
#' @param object Fitted object of class \code{SDTree}.
#' @seealso \code{\link{SDTree}}
plot.SDTree <- function(object){
  # plot function for the spectral deconfounded tree
  data.tree::SetEdgeStyle(object$tree, label = function(x) {x$decision})
  data.tree::SetNodeStyle(object$tree, label = function(x) {x$label})
  plot(object$tree)
}

#### SDForest functions ####
SDForest <- function(formula = NULL, data = NULL, x = NULL, y = NULL, nTree = 100, max_leaves = 500, 
                     cp = 0, min_sample = 3, mtry = NULL, multicore = F, mc.cores = NULL, 
                     Q_type = 'trim', trim_quantile = 0.5, confounding_dim = 0, Q = NULL){

  input_data <- data.handler(formula = formula, data = data, x = x, y = y)
  X <- input_data$X
  Y <- input_data$Y

  # number of observations
  n <- nrow(X)
  # number of covariates
  p <- ncol(X)

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
    mtry <- floor(0.5 * p)
  }

  # bootstrap samples
  ind <- lapply(1:nTree, function(x)sample(1:n, n, replace = T))

  if(multicore){
    if(!is.null(mc.cores)){
      n_cores <- mc.cores
    }
    if(locatexec::is_unix()){
      res <- parallel::mclapply(ind, function(i)SDTree(x = X[i, ], y = Y[i], max_leaves = max_leaves, cp = cp, 
                                            min_sample = min_sample, Q_type = Q_type, 
                                            trim_quantile = trim_quantile, confounding_dim = confounding_dim, mtry = mtry), 
                                            mc.cores = n_cores)
    }else{
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      parallel::clusterExport(cl, c("SDTree", "get_Q", "data.handler", "get_all_splitt", 
                                    "find_s", "evaluate_splitt", "loss", "predict_outsample", 
                                    "traverse_tree", "splitt_names", "leave_names"))
      res <- parallel::clusterApplyLB(cl = cl, i = ind, fun = function(i)SDTree(x = X[i, ], y = Y[i], max_leaves = max_leaves, cp = cp, min_sample = min_sample, 
                  Q_type = Q_type, trim_quantile = trim_quantile, confounding_dim = confounding_dim, mtry = mtry))
      parallel::stopCluster(cl = cl)
    }
  }else{
    res <- pbapply::pblapply(ind, function(i)SDTree(x = X[i, ], y = Y[i], max_leaves = max_leaves, cp = cp, 
                                              min_sample = min_sample, Q_type = Q_type, 
                                              trim_quantile = trim_quantile, confounding_dim = confounding_dim, mtry = mtry))
  }

  # ensemble predictions for each observation
  # but only with the trees that did not contain the observation in the training set
  oob_ind <- lapply(1:n, function(i) which(unlist(lapply(lapply(ind, 
                         function(train)c(1:n)[-train]), 
                         function(x) any(x == i)))))

  oob_predictions <- unlist(lapply(1:n, function(i){
    if(length(oob_ind[[i]]) == 0){
      return(NA)
    }
    predictions <- lapply(oob_ind[[i]], function(model){
      predict_outsample(res[[model]]$tree, X[i, ])
    })
    return(mean(unlist(predictions)))
  }))

  oob_SDloss <- loss(Q %*% Y, Q %*% oob_predictions)
  oob_loss <- loss(Y, oob_predictions)

  # predict with all trees
  pred <- do.call(cbind, lapply(res, function(x){predict_outsample(x$tree, X)}))
  
  # use mean over trees as final prediction
  f_X_hat <- rowMeans(pred)

  # variable importance
  var_imp <- sapply(res, function(x){x$var_importance})
  if(p > 1){
    var_imp <- rowMeans(var_imp)
  }else {
    var_imp <- mean(var_imp)
  }

  output <- list(predictions = f_X_hat, forest = res, var_names = colnames(data.frame(X)), 
                 oob_loss = oob_loss, oob_SDloss = oob_SDloss, var_importance = var_imp, 
                 oob_ind = oob_ind, X = X, Y = Y, Q = Q)
  class(output) <- 'SDForest'
  return(output)
}

predict.SDForest <- function(object, newdata){
  # predict function for the spectral deconfounded random forest
  # using the mean over all trees as the prediction
  # check data type
  if(!is.data.frame(newdata)) stop('newdata must be a data.frame')
  if(!all(object$var_names %in% names(newdata))) stop('newdata must contain all covariates used for training')

  X <- newdata[, object$var_names]
  if(any(is.na(X))) stop('X must not contain missing values')

  pred <- do.call(cbind, lapply(object[[2]], function(x){predict_outsample(x$tree, X)}))
  return(rowMeans(pred))
}



#TODO: permutation importance
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
    if(new_name %in% node$Get('name', filterFun = data.tree::isLeaf)){
        new_name <- paste(new_name, '')
    }
    node$label <- new_name
}

evaluate_splitt <- function(branch, j, s, index, X_branch_j, Y_tilde, Q_temp, n, min_sample){
  # evaluate a split at partition branch on covariate j at the splitpoint s
  # index: index of observations in branch
  # dividing observation in branch
  if (sum(X_branch_j > s) < min_sample | sum(X_branch_j <= s) < min_sample){
    return(list('dloss' = 0, j = j, s = s, branch = branch))
  }

  e_next <- matrix(0, n, 1)
  e_next[index[X_branch_j > s]] <- 1

  u_next_prime <- Q_temp %*% e_next
  u_next_size <- sum(u_next_prime ** 2)

  dloss <- crossprod(u_next_prime, Y_tilde)**2 / u_next_size

  return(list('dloss' = dloss, j = j, s = s, branch = branch))
}

get_all_splitt <- function(branch, X, Y_tilde, Q_temp, n, min_sample, p, E){
  # finds the best splitts for every covariate in branch
  # returns the best splitpoint for every covariate and the resulting loss decrease

  # observations belonging to branch
  index <- which(E[, branch] == 1)

  X_branch <- X[index, ]

  # all possible split points
  s <- find_s(unique(X_branch), min_sample, p)

  res <- lapply(1:p, function(j) {lapply(s[, j], function(x) {
            X_branch_j <- if(p == 1) X_branch else X_branch[, j]
            eval <- evaluate_splitt(branch = branch, j = j, 
              s = x, index = index, X_branch_j = X_branch_j, Y_tilde = Y_tilde, Q_temp = Q_temp, n = n, 
              min_sample = min_sample)
            return(eval)})})

  res <- lapply(res, function(x)do.call(rbind, x))
  res_opt <- lapply(res, function(x)x[which.max(x[, 1]), ])
  res_opt <- do.call(rbind, res_opt)

  return(matrix(unlist(res_opt), ncol = 4, byrow = F))
}

find_s <- function(X, min_sample, p){
  # finds all the reasnable splitting points in a data matrix
  
  if(p == 1){
    X <- matrix(X, ncol = 1)
  }
  n <- nrow(X)

  X_sort <- apply(X, 2, sort, method = 'quick')
  if(min_sample > 1) {
    X_sort <- X_sort[-c(1:(min_sample-1), (n - min_sample + 2):n), ]
  }
  
  if(is.null(dim(X_sort))){
    X_sort <- matrix(X_sort, ncol = 1)
  }
  # find middlepoints between observed x values
  s <- X_sort[-nrow(X_sort), ] + diff(X_sort)/2

  # for runtime reasons
  if(dim(s)[1] > 1000){
    s <- s[seq(1, dim(s)[1], 100), ]
  }else if (dim(s)[1] > 200) {
    s <- s[seq(1, dim(s)[1], 5), ]
  }else if (dim(s)[1] > 100) {
    s <- s[seq(1, dim(s)[1], 2), ]
  }
  
  if(is.null(dim(s))){
    s <- matrix(s, ncol = 1)
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
  if(is.null(dim(X))){
    return(traverse_tree(tree, X))
  }
  return(apply(X, 1, function(x)traverse_tree(tree, x)))
}

loss <- function(Y, f_X){
  # MSE
  return(sum((Y - f_X)^2) / length(Y))
}

pruned_loss <- function(tree, X_val, Y_val, Q_val, t){
  # funtion to prune tree using the minimum loss decrease t
  # and return spectral loss on the validation set
  
  tree_t <- data.tree::Clone(tree)
  
  # prune tree
  data.tree::Prune(tree_t, function(x) x$dloss > t)
  
  # predict on test set
  f_X_hat_val <- predict_outsample(tree_t, X_val)
  
  # return spectral loss
  return(sum((Q_val %*% Y_val - Q_val %*% f_X_hat_val) ** 2) / length(Y_val))
}

varImp <- function(object) {UseMethod("varImp")}

varImp.SDTree <- function(object){
  j_dec <- object$tree$Get(function(x)c(x$j, x$res_dloss), filterFun = function(x)!data.tree::isLeaf(x))
  var_importance <- rep(0, length(object$var_names))
  if(is.null(j_dec)){
    names(var_importance) <- object$var_names
    return(var_importance)
  }
  for(i in 1:ncol(j_dec)){
    var_importance[j_dec[1, i]] <- var_importance[j_dec[1, i]] + j_dec[2, i]
  }
  names(var_importance) <- object$var_names
  return(var_importance)
}

varImp.SDForest <- function(object){
  rowMeans(sapply(object$forest, varImp))
}

prune <- function(object, ...) UseMethod('prune')

prune.SDTree <- function(object, cp){
  object$tree <- data.tree::Clone(object$tree)
  data.tree::Prune(object$tree, function(x) x$cp > cp)
  object$tree$Do(leave_names, filterFun = data.tree::isLeaf)
  object$predictions <- NULL
  object$var_importance <- varImp(object)
  return(object)
}

prune.SDForest <- function(forest, cp, oob = T){
  n <- length(forest$Y)

  pruned_forest <- lapply(forest$forest, function(tree){prune(tree, cp)})
  forest$forest <- pruned_forest

  if(oob){
    oob_predictions <- unlist(lapply(1:n, function(i){
      if(length(forest$oob_ind[[i]]) == 0){
        return(NA)
      }
      predictions <- lapply(forest$oob_ind[[i]], function(model){
        predict_outsample(forest$forest[[model]]$tree, forest$X[i, ])
      })
      return(mean(unlist(predictions)))
    }))
    forest$oob_SDloss <- loss(forest$Q %*% forest$Y, forest$Q %*% oob_predictions)
    forest$oob_loss <- loss(forest$Y, oob_predictions)

    # predict with all trees
    pred <- do.call(cbind, lapply(forest$forest, function(x){predict_outsample(x$tree, forest$X)}))
    
    # use mean over trees as final prediction
    f_X_hat <- rowMeans(pred)
    forest$predictions <- f_X_hat
  }

  # variable importance
  forest$var_importance <- rowMeans(sapply(forest$forest, function(x){x$var_importance}))  

  return(forest)
}

regPath <- function(object, ...) UseMethod('regPath')

regPath.SDTree <- function(object){
  cp_seq <- c(seq(0, 0.1, 0.001), seq(0.1, 0.5, 0.03), seq(0.5, 1, 0.1))
  res <- lapply(cp_seq, function(cp){
    pruned_object <- prune(object, cp)
    return(list(var_importance = pruned_object$var_importance))})

  varImp_path <- t(sapply(res, function(x)x$var_importance))

  paths <- list(cp = cp_seq, varImp_path = varImp_path)
  class(paths) <- 'paths'
  return(paths)
}

regPath.SDForest <- function(object, oob = F, multicore = F, mc.cores = NULL){
  cp_seq <- c(seq(0, 0.1, 0.001), seq(0.1, 0.5, 0.03), seq(0.5, 1, 0.1))
  
  if(multicore){
    if(!is.null(mc.cores)){
      n_cores <- mc.cores
    }
    res <- parallel::mclapply(cp_seq, function(cp){
      pruned_object <- prune(object, cp, oob = oob)
      return(list(var_importance = pruned_object$var_importance, 
                  oob_SDloss = pruned_object$oob_SDloss, 
                  oob_loss = pruned_object$oob_loss))}, 
                  mc.cores = n_cores)
  }else{
    res <- pbapply::pblapply(cp_seq, function(cp){
      pruned_object <- prune(object, cp, oob = oob)
      return(list(var_importance = pruned_object$var_importance, 
                  oob_SDloss = pruned_object$oob_SDloss, 
                  oob_loss = pruned_object$oob_loss))})
  }

  varImp_path <- t(sapply(res, function(x)x$var_importance))
  if(!oob){
    paths <- list(cp = cp_seq, varImp_path = varImp_path)
    class(paths) <- 'paths'
    return(paths)
  }
  
  loss_path <- t(sapply(res, function(x)c(x$oob_SDloss, x$oob_loss)))
  colnames(loss_path) <- c('oob SDE', 'oob MSE')
  paths <- list(cp = cp_seq, varImp_path = varImp_path, loss_path = loss_path, 
                cp_min = cp_seq[which.min(loss_path[, 1])])
  class(paths) <- 'paths'
  return(paths)
}

stabilitySelection <- function(object, ...) UseMethod('stabilitySelection')

stabilitySelection.SDForest <- function(object, multicore = F, mc.cores = NULL){
  cp_seq <- regPath(object$forest[[1]])$cp

  if(multicore){
    if(!is.null(mc.cores)){
      n_cores <- mc.cores
    }
    imp <- parallel::mclapply(object$forest, function(x)regPath(x)$varImp_path > 0, mc.cores = n_cores)
  }else{
    imp <- pbapply::pblapply(object$forest, function(x)regPath(x)$varImp_path > 0)
  }
  imp <- lapply(imp, function(x)matrix(as.numeric(x), ncol = ncol(x)))
  imp <- Reduce('+', imp) / length(object$forest)
  names(imp) <- object$var_names
  paths <- list(cp = cp_seq, varImp_path = imp)
  class(paths) <- 'paths'
  return(paths)
}

plot.paths <- function(object, plotly = F){
  imp_data <- data.frame(object$varImp_path, cp = object$cp)
  imp_data <- tidyr::gather(imp_data, key = 'covariate', value = 'importance', -cp)
  
  gg_path <- ggplot2::ggplot(imp_data, ggplot2::aes(x = cp, y = importance, col = covariate)) +
      ggplot2::geom_line() + 
      ggplot2::theme_bw()

  if(plotly) return(plotly::ggplotly(gg_path))
  return(gg_path)
}

plotOOB <- function(object){
    loss_data <- data.frame(object$loss_path, cp = object$cp)
    gg_sde <- ggplot2::ggplot(loss_data, ggplot2::aes(x = cp, y = oob.SDE)) +
        ggplot2::geom_line() + 
        ggplot2::theme_bw()

    gg_mse <- ggplot2::ggplot(loss_data, ggplot2::aes(x = cp, y = oob.MSE)) +
        ggplot2::geom_line() + 
        ggplot2::theme_bw()
    gridExtra::grid.arrange(gg_sde, gg_mse, ncol = 2)
}

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
      return(list(X = X, Y = as.numeric(Y)))
    }
  }
}

simulate_data_nonlinear <- function(q, p, n, m, eff = NULL){
    #simulate data with confounding and non-linear f_X
    # q: number of confounding covariates in H
    # p: number of covariates in X
    # n: number of observations
    # m: number of covariates with a causal effect on Y

    # complexity of f_X
    complexity <- 5
    # random parameter for fourier basis
    beta <- runif(m * complexity * 2, -1, 1)

    # random confounding covariates H
    H <- matrix(rnorm(n * q, 0, 1), nrow = n)

    # random correlation matrix cov(X, H)
    Gamma <- matrix(rnorm(q * p, 0, 1), nrow = q)

    if(!is.null(eff)){
        non_effected <- p - eff
        if(non_effected <= 0) stop('eff must be smaller than p or NULL')
        
        Gamma[, sample(1:p, non_effected)] <- 0
    }

    # random coefficient vector delta
    delta <- rnorm(q, 0, 1)

    # random error term
    E <- matrix(rnorm(n * p, 0, 1), nrow = n)

    if(q == 0){
        X <- E
    }else{
        X <- H %*% Gamma + E
    }
  
    # random sparse subset of covariates in X
    js <- sample(1:p, m)

    # generate f_X
    f_X <- apply(X, 1, function(x) f_four(x, beta, js))
    
    # generate Y
    Y <- f_X + H %*% delta + rnorm(n, 0, 0.1)
  
    #return data
    return(list(X = X, Y = Y, f_X = f_X, j = js, beta = beta, H = H))
}

f_four <- function(x, beta, js){
    # function to generate f_X
    # x: covariates
    # beta: parameter vector
    # js: relevant covariates

    # number of relevant covariates
    m <- length(js)

    # complexity of f_X
    complexity <- length(beta) / (2 * m)

    # calculate f_X
    do.call(sum, lapply(1:m, function(i) {
        j <- js[i]
        # select beta for covariate j
        beta_ind <- 1:(2*complexity) + (i-1) * 2 * complexity
        # calculate f_X_j
        do.call(sum, lapply(1:complexity, function(k) beta[beta_ind[1 + (k-1) *2]] * sin(k * 0.1 * x[j]) + beta[beta_ind[2 + (k-1) *2]] * cos(k * 0.1 * x[j])))
        }))
}