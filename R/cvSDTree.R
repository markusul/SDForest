#' Cross-validation for the SDTree
#' 
#' Estimates the optimal complexity parameter for the SDTree using cross-validation. 
#' Q is estimated for each training set and validation set separately to ensure independence of the validation set.
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
#' @param mc.cores Number of cores to use for parallel computing, if \code{NULL} all available cores - 1 are used.
#' @param Q_type Type of deconfounding, one of 'trim', 'DDL_trim', 'pca', 'no_deconfounding'. 'trim' corresponds to the Trim transform \insertCite{Cevid2020SpectralModels}{SDForest}, 'DDL_trim' to the implementation of the Doubly debiased lasso \insertCite{Guo2022DoublyConfounding}{SDForest}, 'pca' to the PCA transformation\insertCite{Paul2008PreconditioningProblems}{SDForest} and 'no_deconfounding' to the Identity. See \code{\link{get_Q}}.
#' @param trim_quantile Quantile for Trim transform and DDL Trim transform, only needed for trim and DDL_trim, see \code{\link{get_Q}}.
#' @param q_hat Assumed confounding dimension, only needed for pca, see \code{\link{get_Q}}.
#' @param n_cv Number of folds for cross-validation. It is recommended to not use more than 5 folds if the number of covariates is larger than the number of observations. In this case the spectral transformation could differ to much if the validation data is substantially smaller than the training data.
#' @param cp_seq Sequence of complexity parameters cp to compare using cross-validation, if \code{NULL} a sequence from 0 to 0.6 with stepsize 0.002 is used.
#' @return A list containing
#' \item{cp_min}{The optimal complexity parameter.}
#' \item{cp_table}{A table containing the complexity parameter, the mean and the standard deviation of the loss on the validation sets for the complexity parameters. If multiple complexity parameters result in the same loss, only the one with the largest complexity parameter is shown.}
#' @examples
#' # TODO: add example
#' @seealso \code{\link{SDTree}}
#' @export
cvSDTree <- function(formula = NULL, data = NULL, x = NULL, y = NULL, max_leaves = 50, 
                      min_sample = 5, fast = TRUE, Q_type = 'trim', trim_quantile = 0.5, 
                      q_hat = 0, mc.cores = 1, n_cv = 3, cp_seq = NULL){

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
  Q <- get_Q(X, Q_type, trim_quantile, q_hat)
  
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
    Q_cv <- get_Q(X[cv_ind, ], Q_type, trim_quantile, q_hat)

    X_train <- X[-cv_ind, ]
    Y_train <- Y[-cv_ind]

    X_cv <- X[cv_ind, ]
    Y_cv <- Y[cv_ind]
    
    # estimate tree with the training set
    res <- SDTree(x = X_train, y = Y_train, max_leaves = max_leaves, cp = 0, min_sample = min_sample,
                  Q_type = Q_type, trim_quantile = trim_quantile, q_hat = q_hat)

    # validation performance if we prune with the different ts
    if(mc.cores > 1){
      perf <- parallel::mclapply(t_seq, function(t) pruned_loss(res$tree, X_cv, Y_cv, Q_cv, t), 
                       mc.cores = mc.cores, mc.preschedule = FALSE)
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