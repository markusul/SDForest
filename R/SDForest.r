



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
#' @param Q_type Type of deconfounding, one of 'trim', 'DDL_trim', 'pca', 'no_deconfounding'. 'trim' corresponds to the Trim transform \insertCite{Cevid2020SpectralModels}{SDForest}, 'DDL_trim' to the implementation of the Doubly debiased lasso \insertCite{Guo2022DoublyConfounding}{SDForest}, 'pca' to the PCA transformation\insertCite{Paul2008PreconditioningProblems}{SDForest} and 'no_deconfounding' to the Identity. See \code{\link{get_Q}}.
#' @param trim_quantile Quantile for Trim transform and DDL Trim transform, only needed for trim and DDL_trim, see \code{\link{get_Q}}.
#' @param q_hat Assumed confounding dimension, only needed for pca, see \code{\link{get_Q}}.
#' @param Q Spectral transformation, if \code{NULL} it is internally estimated using \code{\link{get_Q}}.
#' @return Object of class \code{SDTree} containing
#' \item{predictions}{Predictions for the training set.}
#' \item{tree}{The estimated tree of type \code{data.tree} \insertCite{data.tree}{SDForest}. The tree contains the information about all the splits and the resulting estimates.}
#' \item{var_names}{Names of the covariates in the training data.}
#' @examples
#' # TODO: add example
#' @export
SDTree <- function(formula = NULL, data = NULL, x = NULL, y = NULL, max_leaves = NULL, cp = 0.01, min_sample = 5, mtry = NULL, fast = TRUE,
                   Q_type = 'trim', trim_quantile = 0.5, q_hat = 0, Q = NULL, 
                   A = NULL, gamma = 0.5, gpu = FALSE, mem_size = 1e+7, max_candidates = 100){
  ifelse(GPUmatrix::installTorch(), gpu_type <- 'torch', gpu_type <- 'tensorflow')
  input_data <- data.handler(formula = formula, data = data, x = x, y = y)
  X <- input_data$X
  Y <- input_data$Y

  # number of observations
  n <- dim(X)[1]
  # number of covariates
  p <- dim(X)[2]

  if(is.null(max_leaves)) max_leaves <- n


  m <- max_leaves - 1
  mem_size <- mem_size / n
  # check validity of input
  if(n != length(Y)) stop('X and Y must have the same number of observations')
  if(m < 1) stop('max_leaves must be larger than 1')
  if(min_sample < 1) stop('min_sample must be larger than 0')
  if(cp < 0) stop('cp must be at least 0')
  if(!is.null(mtry) && mtry < 1) stop('mtry must be larger than 0')
  if(!is.null(mtry) && mtry > p) stop('mtry must be at most p')
  if(n < 2 * min_sample) stop('n must be at least 2 * min_sample')
  if(max_candidates < 1) stop('max_candidates must be at least 1')

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

  # calculate first estimate
  E <- matrix(1, n, 1)

  E_tilde <- matrix(rowSums(Q))
  if(gpu){
    E_tilde <- gpu.matrix(E_tilde, type = gpu_type)
  }
  
  
  u_start <- E_tilde / sqrt(sum(E_tilde ** 2))
  Q_temp <- Q - u_start %*% (t(u_start) %*% Q)

  Y_tilde <- Q %*% Y
  
  # solve linear model
  if(gpu_type == 'tensorflow'){
    c_hat <- lm.fit(as.matrix(E_tilde), as.matrix(Y_tilde))$coefficients
  }else{
    c_hat <- qr.coef(qr(E_tilde), Y_tilde)
    c_hat <- as.numeric(c_hat)
  }

  #c_hat <- lm.fit(E_tilde, Y_tilde)$coefficients

  loss_start <- as.numeric(sum((Y_tilde - c_hat) ** 2) / n)
  loss_temp <- loss_start

  # initialize tree
  tree <- data.tree::Node$new(name = '1', value = as.numeric(c_hat), 
    dloss = as.numeric(loss_start), cp = 10, n_samples = n)
  
  # memory for optimal splits
  #memory <- lapply(1:(m + 1), function(rep) matrix(0, p, 4))
  memory <- list()
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
      to_small <- sapply(potential_splitts, function(x){sum(E[, x]) < min_sample*2})
      potential_splitts <- potential_splitts[!to_small]
    }

    #iterate over new to estimate splits
    for(branch in potential_splitts){
      E_branch <- E[, branch]
      index <- which(E_branch == 1)
      X_branch <- as.matrix(X[index, ])
    
      s <- find_s(X_branch, max_candidates = max_candidates)
      if(min_sample > 1) {
        s <- s[-c(0:(min_sample - 1), (n - min_sample + 2):(n+1)), ]
      }
      s <- as.matrix(s)

      all_n_splits <- apply(s, 2, function(x) length(unique(x)))
      all_idx <- cumsum(all_n_splits)

      eval <- matrix(-Inf, nrow(s), p)
      done_splits <- 0
      p_top <- 0
      while(p_top < p){
        c_all_idx <- all_idx - done_splits
        p_low <- p_top + 1
        possible <- which(c_all_idx < mem_size)
        p_top <- possible[length(possible)]

        c_n_splits <- sum(all_idx[p_top], -all_idx[p_low-1])
        E_next <- matrix(0, n, c_n_splits)
        for(j in p_low:p_top){
          s_j <- s[, j]
          s_j <- unique(s_j)
          for(i_s in 1:all_n_splits[j]){
            E_next[index[X_branch[, j] > s_j[i_s]], sum(c_all_idx[j-1], i_s)] <- 1
          }
        }
        if(gpu) E_next <- gpu.matrix(E_next, type = gpu_type)

        U_next_prime <- Q_temp %*% E_next
        U_next_size <- colSums(U_next_prime ** 2)
        dloss <- as.numeric(crossprod(U_next_prime, Y_tilde))**2 / U_next_size
        dloss[is.na(dloss)] <- 0
        
        for(m in p_low:p_top){
          eval[1:all_n_splits[m], m] <- dloss[sum(c_all_idx[m-1], 1):c_all_idx[m]]
        }
        done_splits <- done_splits + c_n_splits
      }
      is_opt <- apply(eval, 2, which.max)
      memory[[branch]] <- t(sapply(1:p, function(j) c(eval[is_opt[j], j], j, unique(s[, j])[is_opt[j]], branch)))
    }

    if(i > after_mtry && !is.null(mtry)){
      Losses_dec <- lapply(memory, function(branch){
        branch[sample(1:p, mtry), ]})
      Losses_dec <- do.call(rbind, Losses_dec)
    }else {
       Losses_dec <- do.call(rbind, memory)
    }

    loc <- which.max(Losses_dec[, 1])
    best_branch <- Losses_dec[loc, 4]
    j <- Losses_dec[loc, 2]
    s <- Losses_dec[loc, 3]

    if(Losses_dec[loc, 1] <= 0){
      break
    }

    # divide observations in leave
    index <- which(E[, best_branch] == 1)
    index_n_branches <- index[X[index, j] > s]


    # new indicator matrix
    E <- cbind(E, matrix(0, n, 1))
    E[index_n_branches, best_branch] <- 0
    E[index_n_branches, i+1] <- 1

    E_tilde_branch <- E_tilde[, best_branch]
    suppressWarnings({
    E_tilde[, best_branch] <- Q %*% E[, best_branch]
    })
    E_tilde <- cbind(E_tilde, matrix(E_tilde_branch - E_tilde[, best_branch]))

    if(gpu_type == 'tensorflow'){
      c_hat <- lm.fit(as.matrix(E_tilde), as.matrix(Y_tilde))$coefficients
    }else{
      c_hat <- qr.coef(qr(E_tilde), Y_tilde)
    }

    u_next_prime <- Q_temp %*% E[, i + 1]
    u_next <- u_next_prime / sqrt(sum(u_next_prime ** 2))

    Q_temp <- Q_temp - u_next %*% (t(u_next) %*% Q)

    # check if loss decrease is larger than minimum loss decrease
    # and if linear model could be estimated

    if(sum(is.na(as.matrix(c_hat))) > 0){
      warning('singulaer matrix QE, tree might be to large, consider increasing cp')
      break
    }

    loss_dec <- as.numeric(loss_temp - loss(Y_tilde, E_tilde %*% c_hat))
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
    leave$AddChild(best_branch, value = 0, dloss = loss_dec, cp = loss_dec / loss_start, 
      decision = 'no', n_samples = sum(E[, best_branch] == 1))
    leave$AddChild(i + 1, value = 0, dloss = loss_dec, cp = loss_dec / loss_start, 
      decision = 'yes', n_samples = sum(E[, i + 1] == 1))

    # add estimates to tree leaves
    c_hat <- as.numeric(c_hat)
    for(l in tree$leaves){
      l$value <- c_hat[as.numeric(l$name)]
    }

    # the two new partitions need to be checked for optimal splits in next iteration
    potential_splitts <- c(best_branch, i + 1)

    # a partition with less than min_sample observations or unique samples are not available for further splits
    to_small <- sapply(potential_splitts, function(x){
      new_samples <- nrow(unique(X[as.logical(E[, x]),]))
      if(is.null(new_samples)) new_samples <- 0
      (new_samples < min_sample * 2)
      })
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

  # cp max of all splits after
  tree$Do(function(node) node$cp_max <- max(node$Get('cp')))
  tree$Do(function(node) {
    cp_max <- data.tree::Aggregate(node, 'cp_max', max)
    node$children[[1]]$cp_max <- cp_max
    node$children[[2]]$cp_max <- cp_max
    }, filterFun = data.tree::isNotLeaf
  )

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
cv.SDTree <- function(formula = NULL, data = NULL, x = NULL, y = NULL, max_leaves = 50, 
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
    suppressWarnings({
    res <- SDTree(x = X_train, y = Y_train, max_leaves = max_leaves, cp = 0, min_sample = min_sample,
                  Q_type = Q_type, trim_quantile = trim_quantile, q_hat = q_hat)
    })

    # validation performance if we prune with the different ts
    if(mc.cores > 1){
      perf <- mclapply(t_seq, function(t) pruned_loss(res$tree, X_cv, Y_cv, Q_cv, t), 
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

#' @export 
SDForest <- function(formula = NULL, data = NULL, x = NULL, y = NULL, nTree = 100, 
                     cp = 0, min_sample = 5, mtry = NULL, mc.cores = 1, 
                     Q_type = 'trim', trim_quantile = 0.5, q_hat = 0, Q = NULL, 
                     A = NULL, gamma = 0.5, max_size = NULL, gpu = FALSE, return_data = TRUE, 
                     mem_size = 1e+7, leave_out_ind = NULL, leave_out_envs = NULL, each_trees = NULL, 
                     max_candidates = 100){
                      
  ifelse(GPUmatrix::installTorch(), gpu_type <- 'torch', gpu_type <- 'tensorflow')
  input_data <- data.handler(formula = formula, data = data, x = x, y = y)
  X <- input_data$X
  Y <- input_data$Y

  # number of observations
  n <- nrow(X)
  # number of covariates
  p <- ncol(X)

  if(is.null(max_size)) max_size <- n

  if(n != length(Y)) stop('X and Y must have the same number of observations')
  if(!is.null(mtry) && mtry < 1) stop('mtry must be larger than 0')
  if(!is.null(mtry) && mtry > p) stop('mtry must be at most p')
  if(gpu && (mc.cores > 1)) warning('gpu and multicore cannot be used together, no gpu is not used for tree estimations')
  if(!is.null(leave_out_ind) && !is.null(leave_out_envs)) stop('leave_out_ind and leave_out_envs cannot be used together')
  if(!is.null(leave_out_ind) && any(leave_out_ind > n, leave_out_ind < 1)) stop('leave_out_ind must be smaller than n')
  if(!is.null(leave_out_envs) && length(leave_out_envs) != n) stop('leave_out_envs must have length n')
  if(!is.null(leave_out_envs) && any(is.null(each_trees), each_trees < 1)) stop('each_trees must be larger than 0 if leave_out_envs is used')


  if(!is.null(A)){
    if(is.null(gamma)) stop('gamma must be provided if A is provided')
    if(!is.matrix(A)) stop('A must be a matrix')
    if(nrow(A) != n) stop('A must have n rows')
    W <- get_W(A, gamma, gpu)
  }else {
    W <- diag(n)
    if(gpu) W <- gpu.matrix(W, type = gpu_type)
  }

  # estimate spectral transformation
  if(is.null(Q)){
    Q <- get_Q(as.matrix(W %*% X), Q_type, trim_quantile, q_hat, gpu)
  }else{
    if(!is.matrix(Q)) stop('Q must be a matrix')
    if(any(dim(Q) != n)) stop('Q must have dimension n x n')
  }

  Q <- Q %*% W

  # mtry
  if(is.null(mtry)){
    mtry <- floor(0.5 * p)
    if(mtry < 1) mtry <- 1
  }

  # bootstrap samples
  #TODO: add support for leave_out_envs
  all_ind <- 1:n
  if(!is.null(leave_out_ind)){
    all_ind <- all_ind[-leave_out_ind]
  }
  ind <- lapply(1:nTree, function(x)sample(all_ind, min(length(all_ind), max_size), replace = T))

  if(mc.cores > 1){
    if(locatexec::is_unix()){
      res <- parallel::mclapply(ind, function(i)SDTree(x = X[i, ], y = Y[i], cp = cp, 
                                            min_sample = min_sample, Q_type = Q_type, 
                                            trim_quantile = trim_quantile, q_hat = q_hat, mtry = mtry, 
                                            A = A[i, ], gamma = gamma, mem_size = mem_size), 
                                            max_candidates = max_candidates, 
                                mc.cores = mc.cores)
    }else{
      cl <- parallel::makeCluster(mc.cores)
      doParallel::registerDoParallel(cl)
      parallel::clusterExport(cl, c("SDTree", "get_Q", "data.handler",
                                    "find_s", "loss", "predict_outsample", 
                                    "traverse_tree", "splitt_names", "leave_names"))
      res <- parallel::clusterApplyLB(cl = cl, i = ind, fun = function(i)SDTree(x = X[i, ], y = Y[i], cp = cp, min_sample = min_sample, 
                  Q_type = Q_type, trim_quantile = trim_quantile, q_hat = q_hat, mtry = mtry, 
                  A = A[i, ], gamma = gamma, mem_size = mem_size, max_candidates = max_candidates))
      parallel::stopCluster(cl = cl)
    }
  }else{
    res <- pbapply::pblapply(ind, function(i)SDTree(x = X[i, ], y = Y[i], cp = cp, 
                                              min_sample = min_sample, Q_type = Q_type, 
                                              trim_quantile = trim_quantile, q_hat = q_hat, mtry = mtry, 
                                              A = A[i, ], gamma = gamma, gpu = gpu, mem_size = mem_size, max_candidates = max_candidates))
  }

  # ensemble predictions for each observation
  # but only with the trees that did not contain the observation in the training set
  oob_ind <- lapply(1:n, function(i) which(unlist(lapply(lapply(ind, 
                         function(train)c(1:n)[-train]), 
                         function(x) any(x == i)))))

  oob_predictions <- sapply(1:n, function(i){
    if(length(oob_ind[[i]]) == 0){
      return(NA)
    }
    xi <- X[i, ]
    predictions <- sapply(oob_ind[[i]], function(model){
      predict_outsample(res[[model]]$tree, xi)
    })
    return(mean(predictions))
  })

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
                 oob_ind = oob_ind, oob_predictions = oob_predictions)
  if(return_data){
    output$X <- as.matrix(X)
    output$Y <- as.matrix(Y)
    output$Q <- as.matrix(Q)
  }
  class(output) <- 'SDForest'
  return(output)
}

print.SDForest <- function(object){
  cat("SDForest result\n\n")
  cat("Number of trees: ", length(object$forest), "\n")
  cat("Number of covariates: ", length(object$var_names), "\n")
  if(!is.null(object$oob_loss)){
    cat("OOB loss: ", round(object$oob_loss, 2), "\n")
    cat("OOB spectral loss: ", round(object$oob_SDloss, 2), "\n")
  }
}

predictOOB <- function(object, X = NULL){
  if(is.null(X)){
    X <- object$X
  }
  if(!is.matrix(X)) stop('X must be a matrix and either provided or in the object')

  n <- nrow(X)
  oob_ind <- object$oob_ind

  oob_predictions <- sapply(1:n, function(i){
    if(length(oob_ind[[i]]) == 0){
      return(NA)
    }
    xi <- X[i, ]
    predictions <- sapply(oob_ind[[i]], function(model){
      predict_outsample(object$forest[[model]]$tree, xi)
    })
    return(mean(predictions))
  })
  return(oob_predictions)
}

predict.SDForest <- function(object, newdata){
  # predict function for the spectral deconfounded random forest
  # using the mean over all trees as the prediction
  # check data type
  if(!is.data.frame(newdata)) stop('newdata must be a data.frame')
  if(!all(object$var_names %in% names(newdata))) stop('newdata must contain all covariates used for training')

  X <- as.matrix(newdata[, object$var_names])
  if(any(is.na(X))) stop('X must not contain missing values')

  pred <- do.call(cbind, lapply(object[[2]], function(x){predict_outsample(x$tree, X)}))
  return(rowMeans(pred))
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




mergeForest <- function(fit1, fit2){
  if(any(fit1$var_names != fit2$var_names)) stop('forest must be trained using the same covariates')

  len_1 <- length(fit1$forest)
  len_2 <- length(fit2$forest)
  len_new <- len_1 + len_2

  fit1$forest <- c(fit1$forest, fit2$forest)
  fit1$var_importance <- (fit1$var_importance * len_1 + 
    fit2$var_importance * len_2) / len_new

  if(all(dim(fit1$X) == dim(fit2$X)) && all(fit1$X == fit2$X, fit1$Y == fit2$Y, fit1$Q == fit2$Q, 
    !is.null(fit1$X), !is.null(fit1$Y), !is.null(fit1$Q))){
    fit1$predictions <- (fit1$predictions * length(fit1$forest) + 
      fit2$predictions * length(fit2$forest)) / len_new

    oob_weights_1 <- sapply(fit1$oob_ind, length)
    oob_weights_2 <- sapply(fit2$oob_ind, length)
    oob_weights_new <- oob_weights_1 + oob_weights_2

    oob_predictions_1 <- fit1$oob_predictions * oob_weights_1
    oob_predictions_2 <- fit2$oob_predictions * oob_weights_2

    fit1$oob_predictions <- rowSums(cbind(oob_predictions_1, oob_predictions_2), na.rm = T)
    fit1$oob_predictions <- fit1$oob_predictions / oob_weights_new
    
    fit1$oob_ind <- lapply(1:length(fit1$Y), function(i){
      c(fit1$oob_ind[[i]], fit2$oob_ind[[i]] + len_1)
    })
    fit1$oob_SDloss <- loss(fit1$Q %*% fit1$Y, fit1$Q %*% fit1$oob_predictions)
    fit1$oob_loss <- loss(fit1$Y, fit1$oob_predictions)
  }else{
    warning('forests migth be trained on different data')
    fit1$X <- NULL
    fit1$Y <- NULL
    fit1$Q <- NULL
    fit1$predictions <- NULL
    fit1$oob_prediction <- NULL
    fit1$oob_ind <- NULL
    fit1$oob_loss <- NULL
    fit1$oob_SDloss <- NULL
  }
  return(fit1)
}

toList <- function(object, ...) UseMethod('toList')
fromList <- function(object, ...) UseMethod('fromList')

toList.SDTree <- function(tree){
  tree$tree <- as.list(tree$tree)
  return(tree)
}

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
#TODO: autoForest adding trees until prediction does not change anymore