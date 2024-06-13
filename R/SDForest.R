#' Spectral Deconfounded Random Forest
#' 
#' Estimate regression Random Forest using spectral deconfounding.
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
#' @param nTree Number of trees to grow.
#' @param cp Complexity parameter, minimum loss decrease to split a node. 
#' A split is only performed if the loss decrease is larger than \code{cp * initial_loss}, 
#' where \code{initial_loss} is the loss of the initial estimate using only a stump.
#' @param min_sample Minimum number of observations per leaf. 
#' A split is only performed if both resulting leaves have at least 
#' \code{min_sample} observations.
#' @param mtry Number of randomly selected covariates to consider for a split, 
#' if \code{NULL} half of the covariates are available for each split. 
#' \eqn{\text{mtry} = \lfloor \frac{p}{2} \rfloor}
#' @param mc.cores Number of cores to use for parallel processing,
#' if \code{mc.cores > 1} the trees are estimated in parallel.
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
#' @param max_size Maximum number of observations used for a bootstrap sample.
#' If \code{NULL} n samples with replacement are drawn.
#' @param gpu If \code{TRUE}, the calculations are performed on the GPU. 
#' If it is properly set up.
#' @param return_data If \code{TRUE}, the training data is returned in the output.
#' This is needed for \code{\link{prune.SDForest}}, \code{\link{regPath.SDForest}}, 
#' and for \code{\link{mergeForest}}.
#' @param mem_size Amount of split candidates that can be evaluated at once.
#' This is a trade-off between memory and speed can be decreased if either
#' the memory is not sufficient or the gpu is to small.
#' @param leave_out_ind Indices of observations that should not be used for training.
#' @param envs Vector of environments which can be used for stratified tree fitting.
#' @param leave_envs_out_trees Number of trees that should be estimated while leaving
#' one of the environments out. Results in number of environments times number of trees.
#' NOT SUPPORTED YET
#' @param envs_trees Number of trees that should be estimated for each environment.
#' Results in number of environments times number of trees.
#' NOT SUPPORTED YET
#' @param max_candidates Maximum number of split points that are 
#' proposed at each node for each covariate.
#' @return Object of class \code{SDForest} containing:
#' \item{predictions}{Vector of predictions for each observation.}
#' \item{forest}{List of SDTree objects.}
#' \item{var_names}{Names of the covariates.}
#' \item{oob_loss}{Out-of-bag loss. MSE}
#' \item{oob_SDloss}{Out-of-bag loss using the spectral transformation.}
#' \item{var_importance}{Variable importance.}
#' \item{oob_ind}{List of indices of trees that did not contain the observation in the training set.}
#' \item{oob_predictions}{Out-of-bag predictions.}
#' If \code{return_data} is \code{TRUE} the following are also returned:
#' \item{X}{Matrix of covariates.}
#' \item{Y}{Vector of responses.}
#' \item{Q}{Spectral transformation matrix.}
#' @seealso \code{\link{get_Q}}, \code{\link{get_W}}, \code{\link{SDTree}}, 
#' \code{\link{simulate_data_nonlinear}}, \code{\link{regPath}}, 
#' \code{\link{stabilitySelection}}, \code{\link{prune}}, \code{\link{partDependence}}
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
#' # comparison to classical random forest
#' fit_ranger <- ranger::ranger(Y ~ ., train_data, importance = 'impurity')
#' 
#' fit <- SDForest(x = X, y = Y, nTree = 10, Q_type = 'pca', q_hat = 2)
#' fit <- SDForest(Y ~ ., train_data)
#' fit
#' 
#' # comparison of variable importance
#' imp_ranger <- fit_ranger$variable.importance
#' imp_sdf <- fit$var_importance
#' imp_col <- rep('black', length(imp_ranger))
#' imp_col[sim_data$j] <- 'red'
#' 
#' plot(imp_ranger, imp_sdf, col = imp_col, pch = 20,
#'      xlab = 'ranger', ylab = 'SDForest', 
#'      main = 'Variable Importance')
#' 
#' # check regularization path of variable importance
#' path <- regPath(fit)
#' # out of bag error for different regularization
#' plotOOB(path)
#' plot(path)
#' 
#' # detection of causal parent using stability selection
#' stablePath <- stabilitySelection(fit)
#' plot(stablePath)
#' 
#' # pruning of forest according to optimal out-of-bag performance
#' fit <- prune(fit, cp = path$cp_min)
#' 
#' # partial functional dependence of y on the first causal parent
#' dep <- partDependence(fit, sim_data$j[1])
#' plot(dep, n_examples = 100)
#' @export
SDForest <- function(formula = NULL, data = NULL, x = NULL, y = NULL, nTree = 100, 
                     cp = 0, min_sample = 5, mtry = NULL, mc.cores = 1, 
                     Q_type = 'trim', trim_quantile = 0.5, q_hat = 0, Q = NULL, 
                     A = NULL, gamma = 7, max_size = NULL, gpu = FALSE, 
                     return_data = TRUE, mem_size = 1e+7, leave_out_ind = NULL, 
                     envs = NULL, leave_envs_out_trees = NULL, envs_trees = NULL, 
                     max_candidates = 100){
  if(gpu) ifelse(GPUmatrix::installTorch(), 
                 gpu_type <- 'torch', 
                 gpu_type <- 'tensorflow')
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
  if(gpu && (mc.cores > 1)) 
    warning('gpu and multicore cannot be used together, 
            no gpu is not used for tree estimations')
  #if(!is.null(leave_out_ind) && !is.null(leave_out_envs)) 
  #  stop('leave_out_ind and leave_out_envs cannot be used together')
  #if(!is.null(leave_out_ind) && any(leave_out_ind > n, leave_out_ind < 1)) 
  #  stop('leave_out_ind must be smaller than n')
  #if(!is.null(leave_out_envs) && length(leave_out_envs) != n) 
  #  stop('leave_out_envs must have length n')
  #if(!is.null(leave_out_envs) && any(is.null(each_trees), each_trees < 1)) 
  #  stop('each_trees must be larger than 0 if leave_out_envs is used')

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
  #TODO: add support for stratified fitting
  all_ind <- 1:n
  if(!is.null(leave_out_ind)){
    all_ind <- all_ind[-leave_out_ind]
  }
  ind <- lapply(1:nTree, function(x)
    sample(all_ind, min(length(all_ind), max_size), replace = TRUE))
  
  if(mc.cores > 1){
    if(locatexec::is_unix()){
      print('mclapply')
      res <- parallel::mclapply(ind, function(i) {
        SDTree(x = X[i, ], y = Y[i], cp = cp, min_sample = min_sample, 
               Q_type = Q_type, trim_quantile = trim_quantile, q_hat = q_hat, 
               mtry = mtry, A = A[i, ], gamma = gamma, mem_size = mem_size, 
               max_candidates = max_candidates)
        }, 
        mc.cores = mc.cores)
    }else{
      print('makeCluster')
      cl <- parallel::makeCluster(mc.cores)
      doParallel::registerDoParallel(cl)
      parallel::clusterExport(cl, c("SDTree", "get_Q", "data.handler",
                                    "find_s", "loss", "predict_outsample", 
                                    "traverse_tree", "splitt_names", 
                                    "leave_names", 'X', 'Y', 'A'))
      res <- parallel::clusterApplyLB(cl = cl, ind, fun = function(i)
        SDTree(x = X[i, ], y = Y[i], cp = cp, min_sample = min_sample, 
               Q_type = Q_type, trim_quantile = trim_quantile, q_hat = q_hat, 
               mtry = mtry, A = A[i, ], gamma = gamma, mem_size = mem_size, 
               max_candidates = max_candidates))
      parallel::stopCluster(cl = cl)
    }
  }else{
    res <- pbapply::pblapply(ind, function(i)
      SDTree(x = X[i, ], y = Y[i], cp = cp, min_sample = min_sample, 
             Q_type = Q_type, trim_quantile = trim_quantile, q_hat = q_hat, 
             mtry = mtry, A = A[i, ], gamma = gamma, gpu = gpu, 
             mem_size = mem_size, max_candidates = max_candidates))
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

  output <- list(predictions = f_X_hat, 
                 forest = res, 
                 var_names = colnames(data.frame(X)), 
                 oob_loss = oob_loss, 
                 oob_SDloss = oob_SDloss, 
                 var_importance = var_imp, 
                 oob_ind = oob_ind, 
                 oob_predictions = oob_predictions)
  
  if(return_data){
    output$X <- as.matrix(X)
    output$Y <- as.matrix(Y)
    output$Q <- as.matrix(Q)
  }
  class(output) <- 'SDForest'
  
  output
}
