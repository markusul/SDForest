#' Cross-validation for the SDTree
#' 
#' Estimates the optimal complexity parameter for the SDTree using cross-validation. 
#' The transformations are estimated for each training set and validation set 
#' separately to ensure independence of the validation set.
#' @author Markus Ulmer
#' @references
#'  \insertAllCited{}
#' @param formula Object of class \code{formula} or describing the model to fit 
#' of the form \code{y ~ x1 + x2 + ...} where \code{y} is a numeric response and 
#' \code{x1, x2, ...} are vectors of covariates. Interactions are not supported.
#' @param data Training data of class \code{data.frame} containing the variables in the model.
#' @param x Predictor data, alternative to \code{formula} and \code{data}.
#' @param y Response vector, alternative to \code{formula} and \code{data}.
#' @param max_leaves Maximum number of leaves for the grown tree.
#' @param cp Complexity parameter, minimum loss decrease to split a node. 
#' A split is only performed if the loss decrease is larger than \code{cp * initial_loss}, 
#' where \code{initial_loss} is the loss of the initial estimate using only a stump.
#' @param min_sample Minimum number of observations per leaf. 
#' A split is only performed if both resulting leaves have at least 
#' \code{min_sample} observations.
#' @param mtry Number of randomly selected covariates to consider for a split, 
#' if \code{NULL} all covariates are available for each split.
#' @param fast If \code{TRUE}, only the optimal splitts in the new leaves are 
#' evaluated and the previously optimal splitts and their potential loss-decrease are reused. 
#' If \code{FALSE} all possible splitts in all the leaves are reevaluated after every split.
#' @param Q_type Type of deconfounding, one of 'trim', 'pca', 'no_deconfounding'. 
#' 'trim' corresponds to the Trim transform \insertCite{Cevid2020SpectralModels}{SDForest} 
#' as implemented in the Doubly debiased lasso \insertCite{Guo2022DoublyConfounding}{SDForest}, 
#' 'pca' to the PCA transformation\insertCite{Paul2008PreconditioningProblems}{SDForest}. 
#' See \code{\link{get_Q}}.
#' @param trim_quantile Quantile for Trim transform, 
#' only needed for trim and DDL_trim, see \code{\link{get_Q}}.
#' @param q_hat Assumed confounding dimension, only needed for pca, 
#' see \code{\link{get_Q}}.
#' @param Q Spectral transformation, if \code{NULL} 
#' it is internally estimated using \code{\link{get_Q}}.
#' @param A Numerical Anchor of class \code{matrix}. See \code{\link{get_W}}.
#' @param gamma Strength of distributional robustness, \eqn{\gamma \in [0, \infty]}. 
#' See \code{\link{get_W}}.
#' @param gpu If \code{TRUE}, the calculations are performed on the GPU. 
#' If it is properly set up.
#' @param mem_size Amount of split candidates that can be evaluated at once.
#' This is a trade-off between memory and speed can be decreased if either
#' the memory is not sufficient or the gpu is to small.
#' @param max_candidates Maximum number of split points that are 
#' proposed at each node for each covariate.
#' @param n_cv Number of folds for cross-validation. 
#' It is recommended to not use more than 5 folds if the number of covariates 
#' is larger than the number of observations. In this case the spectral 
#' transformation could differ to much if the validation data is 
#' substantially smaller than the training data.
#' @param cp_seq Sequence of complexity parameters cp to compare using cross-validation, 
#' if \code{NULL} a sequence from 0 to 0.6 with stepsize 0.002 is used.
#' @param mc.cores Number of cores to use for parallel computation.
#' @return A list containing
#' \item{cp_min}{The optimal complexity parameter.}
#' \item{cp_table}{A table containing the complexity parameter, 
#' the mean and the standard deviation of the loss on the validation sets for the 
#' complexity parameters. If multiple complexity parameters result in the same loss, 
#' only the one with the largest complexity parameter is shown.}
#' @examples
#' set.seed(1)
#' n <- 50
#' X <- matrix(rnorm(n * 5), nrow = n)
#' y <- sign(X[, 1]) * 3 + rnorm(n, 0, 5)
#' cp <- cvSDTree(x = X, y = y, Q_type = 'no_deconfounding')
#' cp
#' @seealso \code{\link{SDTree}} \code{\link{prune.SDTree}} \code{\link{regPath.SDTree}}
#' @export
cvSDTree <- function(formula = NULL, data = NULL, x = NULL, y = NULL, 
                     max_leaves = NULL, cp = 0.01, min_sample = 5, mtry = NULL, 
                     fast = TRUE, Q_type = 'trim', trim_quantile = 0.5, q_hat = 0, 
                     Q = NULL, A = NULL, gamma = 0.5, gpu = FALSE, mem_size = 1e+7, 
                     max_candidates = 100, n_cv = 3, cp_seq = NULL, mc.cores = 1){
  ifelse(GPUmatrix::installTorch(), gpu_type <- 'torch', gpu_type <- 'tensorflow')
  input_data <- data.handler(formula = formula, data = data, x = x, y = y)
  X <- input_data$X
  Y <- input_data$Y

  p <- ncol(X)
  n <- nrow(X)

  if(is.null(max_leaves)) max_leaves <- n
  m <- max_leaves - 1
  # check validity of input
  if(n != length(Y)) stop('X and Y must have the same number of observations')
  if(m < 1) stop('max_leaves must be larger than 1')
  if(min_sample < 1) stop('min_sample must be larger than 0')
  if(!is.null(cp_seq) && (any(cp_seq < 0) | any(cp_seq > 1))) 
    stop('cp.seq must be between 0 and 1')
  if(n_cv < 2 | n_cv > n) stop('n_cv must be at least 2 and smaller than n')
  if(n_cv > 5 & n < p) 
    warning('if n < p and to many folds are used, Q_validation might differ 
            to much from Q_trainig, consider using less folds')

  # estimate spectral transformation
  if(!is.null(A)){
    if(is.null(gamma)) stop('gamma must be provided if A is provided')
    if(is.vector(A)) A <- matrix(A)
    if(!is.matrix(A)) stop('A must be a matrix')
    if(nrow(A) != n) stop('A must have n rows')
    W <- get_W(A, gamma, gpu)
  }else {
    
    W <- diag(n)
    if(gpu) W <- gpu.matrix(W, type = gpu_type)
  }

  if(is.null(Q)){
    Q <- get_Q(as.matrix(W %*% X), Q_type, trim_quantile, q_hat, gpu)
  }else{
    if(!is.matrix(Q)) stop('Q must be a matrix')
    if(any(dim(Q) != n)) stop('Q must have dimension n x n')
  }
  Q <- Q %*% W
  
  # estimating initial loss with only a stump
  # to map optimal minimal loss decrease to a cp value
  E <- matrix(1, n, 1)
  E_tilde <- matrix(rowSums(Q))
  if(gpu){
    E_tilde <- gpu.matrix(E_tilde, type = gpu_type)
  }
  Y_tilde <- Q %*% Y

  # solve linear model
  if(gpu_type == 'tensorflow'){
    c_hat <- lm.fit(as.matrix(E_tilde), as.matrix(Y_tilde))$coefficients
  }else{
    c_hat <- qr.coef(qr(E_tilde), Y_tilde)
    c_hat <- as.numeric(c_hat)
  }
  loss_start <- sum((Y_tilde - c_hat) ** 2) / n

  # validation set size
  len_test <- floor(n / n_cv)
  # validation sets
  test_ind <- lapply(1:n_cv, function(x) 1:len_test + (x - 1) * len_test)

  # sequence of minimum loss decrease to compare with cv
  if(is.null(cp_seq)){
    cp_seq <- seq(0, 0.6, 0.002)
  }
  t_seq <- cp_seq * loss_start

  # estimate performance for every validation set
  perf <- lapply(test_ind, function(cv_ind){
    # estimate spectral transformation
    if(!is.null(A)){
      if(is.null(gamma)) stop('gamma must be provided if A is provided')
      if(is.vector(A)) A <- matrix(A)
      if(!is.matrix(A)) stop('A must be a matrix')
      if(nrow(A) != n) stop('A must have n rows')
      W_cv <- get_W(A[cv_ind, ], gamma, gpu)
    }else {
      W_cv <- diag(length(cv_ind))
      if(gpu) W_cv <- gpu.matrix(W_cv, type = gpu_type)
    }
    Q_cv <- get_Q(as.matrix(W_cv %*% X[cv_ind, ]), 
                  Q_type, trim_quantile, q_hat, gpu)
    Q_cv <- Q_cv %*% W_cv

    X_train <- X[-cv_ind, ]
    Y_train <- Y[-cv_ind]
    A_train <- A[-cv_ind, ]

    X_cv <- X[cv_ind, ]
    Y_cv <- Y[cv_ind]

    # estimate tree with the training set
    res <- SDTree(x = X_train, y = Y_train, max_leaves = max_leaves, cp = 0, 
                  min_sample = min_sample, Q_type = Q_type, 
                  trim_quantile = trim_quantile, q_hat = q_hat, mtry = mtry, 
                  A = A_train, gamma = gamma, gpu = gpu, mem_size = mem_size, 
                  max_candidates = max_candidates)

    # validation performance if we prune with the different ts
    if(mc.cores > 1){
      perf <- parallel::mclapply(t_seq, function(t) 
        pruned_loss(res$tree, X_cv, Y_cv, Q_cv, t), 
        mc.cores = mc.cores, mc.preschedule = FALSE)
    }else{
      perf <- lapply(t_seq, function(t) 
        pruned_loss(res$tree, X_cv, Y_cv, Q_cv, t))
    }
    
    return(perf)
  })
  
  # collect performance for different min loss decreases
  perf <- matrix(unlist(perf), ncol = n_cv, byrow = FALSE)
  
  cp_table <- matrix(c(t_seq / loss_start, apply(perf, 1, mean), 
                       apply(perf, 1, sd)), ncol = 3, byrow = FALSE)
  colnames(cp_table) <- c('cp', 'SDLoss mean', 'SDLoss sd')

  loss_unique <- unique(cp_table[, 2])
  cp_table <- lapply(loss_unique, function(loss){
    idx <- which(cp_table[, 2] == loss)
    cp_table[idx[length(idx)], ]
  })
  cp_table <- do.call(rbind, cp_table)

  cp_min <- cp_table[which.min(cp_table[, 2]), 1]

  list(cp_min = cp_min, cp_table = cp_table)
}