
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


