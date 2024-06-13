#' Spectral Deconfounded Tree
#' 
#' Estimates a regression tree using spectral deconfounding. 
#' # TODO: add more details
#' @references
#'  \insertAllCited{}
#' @author Markus Ulmer
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
#' @return Object of class \code{SDTree} containing
#' \item{predictions}{Predictions for the training set.}
#' \item{tree}{The estimated tree of class \code{Node} from \insertCite{Glur2023Data.tree:Structure}{SDForest}. 
#' The tree contains the information about all the splits and the resulting estimates.}
#' \item{var_names}{Names of the covariates in the training data.}
#' \item{var_importance}{Variable importance of the covariates. see \code{\link{varImp.SDTree}}}
#' @seealso \code{\link{simulate_data_nonlinear}}, \code{\link{regPath.SDTree}}, 
#' \code{\link{prune.SDTree}}, \code{\link{partDependence}}
#' @examples
#' set.seed(42)
#' # simulation of confounded data
#' sim_data <- simulate_data_nonlinear(q = 2, p = 150, n = 100, m = 2)
#' X <- sim_data$X
#' Y <- sim_data$Y
#' train_data <- data.frame(X, Y)
#' # causal parents of y
#' sim_data$j
#' 
#' tree_plain_cv <- cvSDTree(Y ~ ., train_data, Q_type = "no_deconfounding")
#' tree_plain <- SDTree(Y ~ ., train_data, Q_type = "no_deconfounding", cp = 0)
#' 
#' tree_causal_cv <- cvSDTree(Y ~ ., train_data)
#' tree_causal <- SDTree(y = Y, x = X, cp = 0)
#' 
#' # check regularization path of variable importance
#' path <- regPath(tree_causal)
#' plot(path)
#' 
#' tree_plain <- prune(tree_plain, cp = tree_plain_cv$cp_min)
#' tree_causal <- prune(tree_causal, cp = tree_causal_cv$cp_min)
#' plot(tree_causal)
#' 
#' 
#' plot(tree_plain)
#' @export
SDTree <- function(formula = NULL, data = NULL, x = NULL, y = NULL, max_leaves = NULL, 
                   cp = 0.01, min_sample = 5, mtry = NULL, fast = TRUE,
                   Q_type = 'trim', trim_quantile = 0.5, q_hat = 0, Q = NULL, 
                   A = NULL, gamma = 0.5, gpu = FALSE, mem_size = 1e+7, max_candidates = 100){
  ifelse(GPUmatrix::installTorch(), gpu_type <- 'torch', gpu_type <- 'tensorflow')
  input_data <- data.handler(formula = formula, data = data, x = x, y = y)
  X <- input_data$X
  Y <- input_data$Y

  # number of observations
  n <- nrow(X)
  # number of covariates
  p <- ncol(X)

  if(is.null(max_leaves)) max_leaves <- n

  max_leaves <- max_leaves - 1

  mem_size <- mem_size / n
  # check validity of input
  if(n != length(Y)) stop('X and Y must have the same number of observations')
  if(max_leaves < 0) stop('max_leaves must be larger than 1')
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

  loss_start <- as.numeric(sum((Y_tilde - c_hat) ** 2) / n)
  loss_temp <- loss_start

  # initialize tree
  tree <- data.tree::Node$new(name = '1', value = as.numeric(c_hat), 
                              dloss = as.numeric(loss_start), 
                              cp = 10, n_samples = n)

  # memory for optimal splits
  memory <- list()
  potential_splitts <- 1

  # variable importance
  var_imp <- rep(0, p)
  names(var_imp) <- colnames(X)

  after_mtry <- 0

  print('a')
  for(i in 1:max_leaves){
    # iterate over all possible splits every time
    # for slow but slightly better solution
    if(!fast){
      potential_splitts <- 1:i
      to_small <- sapply(potential_splitts, 
                         function(x){sum(E[, x]) < min_sample*2})
      potential_splitts <- potential_splitts[!to_small]
    }

    #iterate over new to estimate splits
    for(branch in potential_splitts){
      E_branch <- E[, branch]
      index <- which(E_branch == 1)
      X_branch <- as.matrix(X[index, ])

      s <- find_s(X_branch, max_candidates = max_candidates)
      n_splits <- nrow(s)
      if(min_sample > 1) {
        s <- s[-c(0:(min_sample - 1), (n_splits - min_sample + 2):(n_splits+1)), ]
      }
      s <- matrix(s, ncol = p)

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
      memory[[branch]] <- t(sapply(1:p, function(j) 
        c(eval[is_opt[j], j], j, unique(s[, j])[is_opt[j]], branch)))
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
    leave$AddChild(best_branch, value = 0, dloss = loss_dec, 
                   cp = loss_dec / loss_start, decision = 'no', 
                   n_samples = sum(E[, best_branch] == 1))
    leave$AddChild(i + 1, value = 0, dloss = loss_dec, 
                   cp = loss_dec / loss_start, decision = 'yes', 
                   n_samples = sum(E[, i + 1] == 1))

    # add estimates to tree leaves
    c_hat <- as.numeric(c_hat)
    for(l in tree$leaves){
      l$value <- c_hat[as.numeric(l$name)]
    }

    # the two new partitions need to be checked for optimal splits in next iteration
    potential_splitts <- c(best_branch, i + 1)

    # a partition with less than min_sample observations or unique samples 
    # are not available for further splits
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
    print(i)
  }

  if(i == max_leaves){
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

  res <- list(predictions = f_X_hat, tree = tree, 
              var_names = var_names, var_importance = var_imp)
  class(res) <- 'SDTree'
  
  res
}