library(parallel)
library(data.tree)
library(RcppEigen)
library(data.table)
library(Rdpack)
#' @importFrom Rdpack reprompt

# number of available cores
n_cores <- detectCores()
# if there are less than 24 cores, it will be a local machine leave two cores for other tasks
if(n_cores > 1 && n_cores <= 24){
    n_cores <- n_cores - 1
}

#' Estimation of spectral transformation
#' 
#' Estimates Q
#' \insertCite{Breiman1996BaggingPredictors}{SDForest}
#' @references
#'   \insertRef{Breiman1996BaggingPredictors}{SDForest}
#' @param X numerical covariates of class \code{matrix}
#' @param type type of deconfounding, one of 'trim', 'DDL_trim', 'pca', 'no_deconfounding'
#' @param q_hat assumed confounding dimension, only needed for pca
#' @return Q of class \code{matrix}, the spectral transformation
#' @examples
#' X <- matrix(rnorm(100 * 20), nrow = 100)
#' get_Q(X, 'trim')
#' get_Q(X, 'DDL_trim')
#' get_Q(X, 'pca', q_hat = 5)
#' get_Q(X, 'no_deconfounding')
#' @export
get_Q <- function(X, type, q_hat = 0){
    # X: covariates
    # type: type of deconfounding
    modes <- c('trim' = 1, 'DDL_trim' = 2, 'pca' = 3, 'no_deconfounding' = 4)
    if(!(type %in% names(modes))) stop(paste("type must be one of:", paste(names(modes), collapse = ', ')))

    # number of observations
    n <- dim(X)[1]

    # calculate deconfounding matrix
    sv <- svd(X)
    tau <- median(sv$d)
    D_tilde <- unlist(lapply(sv$d, FUN = function(x)min(x, tau))) / sv$d

    Q <- switch(modes[type], sv$u %*% diag(D_tilde) %*% t(sv$u), # trim
                            diag(n) - sv$u %*% diag(1 - D_tilde) %*% t(sv$u), # DDL_trim
                            { # pca
                                d_pca <- sv$d
                                if(q_hat <= 0) stop("the assumed confounding dimension q_hat must be larger than zero")
                                d_pca[1:q_hat] <- 0
                                sv$u %*% diag(d_pca) %*% t(sv$u)
                            },
                            diag(n)) # no_deconfounding
    return(Q)
}

cv.SDTree <- function(X, Y, m = 100, cp = 0.00001, min_sample = 5, deconfounding = T, multicore = F, n_cv = 3, Q_type = 1){
  # cross-validation to selecet optimal cp
  # m: maximum number of splits, maximum m + 1 leaves
  # cp: cost-complexity parameter
  # min_sample: min number of samples per leaf to split
  # n_cv: number of cross-validation sets
  
  p <- ncol(X)
  n <- length(Y)

  # estimate spectral transformation
  Q <- get_Q(X, Q_type)
  
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
  t_seq <- seq(0, loss_start, 0.1)

  # estimate performance for every validation set
  perf <- lapply(test_ind, function(cv_ind){
    # calculate Trim transform
    Q_cv <- get_Q(X[cv_ind, ], Q_type)

    # deconfound X and Y
    X_train <- X[-cv_ind, ]
    Y_train <- Y[-cv_ind]

    X_cv <- X[cv_ind, ]
    Y_cv <- Y[cv_ind]
    
    # estimate tree with the training set
    res <- SDTree(X_train, Y_train, m = m, cp = cp, min_sample = min_sample, 
                  deconfounding = deconfounding, mtry = FALSE, multicore = multicore, 
                  pruning = F, fast = T, Q_type = Q_type)
    
    # validation performance if we prune with the different ts
    if(multicore){
      perf <- mclapply(t_seq, function(t) val_loss(res$tree, X_cv, Y_cv, Q_cv, t), mc.cores = n_cores)
    }else{
      perf <- lapply(t_seq, function(t) val_loss(res$tree, X_cv, Y_cv, Q_cv, t))
    }
    
    return(perf)
  })
  
  # collect performance for different min loss decreases
  perf <- matrix(unlist(perf), ncol = n_cv, byrow = FALSE)
  # use mean over all cv-sets
  perf <- apply(perf, 1, mean)
  # select minimum
  t.min <- t_seq[[which.min(perf)]]
  
  # return cp.min
  return(t.min / loss_start)
}

val_loss <- function(tree, X_val, Y_val, Q_val, t){
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

SDTree <- function(X, Y, m = 50, cp = 0.00001, min_sample = 5, deconfounding = T, mtry = F, fast = T, multicore = F, pruning = F, val_ratio = 0.3, Q_type = 1){
  # X: covariates
  # Y: outcome
  # m: max number of splits
  # min_sample: min number of samples per leaf to split
  # cp: complexity parameter
  # mtry: if a number, only a random set of size mtry of the covariates are usable for a split
  # pruning: if True the val ratio is used as a validation set and cps down to cp are compared and pruned using the test set
  
  if(pruning){
    # dividing data into train and validation set
    ind <- sample(1:length(Y), ceiling((length(Y)*val_ratio)), replace = F)

    X_val <- X[ind, ]
    Y_val <- Y[ind]

    X <- X[-ind, ]
    Y <- Y[-ind]

    if(deconfounding){
        Q_val <- get_Q(X_val, Q_type)
    }else{
        Q_val <- diag(nrow(X_val))
    }
  }

  # number of observations
  n <- dim(X)[1]
  # number of covariates
  p <- dim(X)[2]

  # calculate deconfounding matrix
  if(deconfounding){
    # calculate trim transform
    Q <- get_Q(X, Q_type)
  }else{
    # no deconfounding
    Q <- diag(n)
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
  if(mtry){
    len_p <- mtry
  }else {
    len_p <- p
  }

  for(i in 1:m){
    if(mtry){
      # sample available covariates
      available_j <- sample(1:p, mtry)
    }else{
      available_j <- 1:p
    }
    
    # iterate over all possible splits every time 
    # for slow but slightly better solution
    if(!fast){
      potential_splitts <- 1:i
      to_small <- unlist(lapply(potential_splitts, function(x){sum(E[, x] == 1) < min_sample}))
      potential_splitts <- potential_splitts[!to_small]
    }
    
    #iterate over new to estimate splits  
    for(branch in potential_splitts){
      best_splitts <- get_all_splitt(branch, X, Y_tilde, Q, n, i+1, E, multicore, min_sample)
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
    index_branch <- index[X[index, j] <= s]
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
    if((Losses_dec[loc] <= cp * loss_start) | (sum(is.na(c_hat)) > 0)){
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
    leave$AddChild(best_branch, value = 0, dloss = Losses_dec[loc])
    leave$AddChild(i + 1, value = 0, dloss = Losses_dec[loc])

    # add estimates to tree leaves
    for(l in tree$leaves){
      l$value <- c_hat[as.numeric(l$name)]
    }
  
    # new temporary loss
    loss_temp <- loss(Y_tilde, E_tilde %*% c_hat)
    
    # the two new partitions need to be checked for optimal splits in next iteration
    potential_splitts <- c(best_branch, i + 1)
    
    # a partition with less than min_sample observations are not available for further splits
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
    print('warning, maximum number of iterations was reached!')
  }
  
  # if pruning is selected, different cp values are tested with the validation and 
  # tree is automatically pruned
  if(pruning){
    # get all the loss decreases for every split in the tree sorted
    s_dloss <- sort(unique(tree$Get('dloss', traversal = 'post-order')))
    # get loss decreases between the ones, that are present in the tree
    # those are the relevant stopping criterias that have to be compared
    t_loss <- s_dloss[-length(s_dloss)] + diff(s_dloss) / 2
    
    # compare all the loss decreases as stopping crition and choose the one with minimum spectral loss
    t <- t_loss[which.min(unlist(lapply(t_loss, function(t)val_loss(tree, X_val, Y_val, Q_val, t))))]
    
    # prune tree with resulting optimal minimum loss decrease
    Prune(tree, function(x) x$dloss > t)
    
    # calculate optimal cp value
    cp_min <- t / loss_start
  }else{
    cp_min <- NA
  }
  
  # predict the test set
  f_X_hat <- predict_outsample(tree, X)

  res <- list(f_X_hat = f_X_hat, tree = tree, cp_min = cp_min)
  class(res) <- 'SDTree'
  return(res)
}

evaluate_splitt <- function(branch, j, s, index, X, Y_tilde, Q = diag(length(Y)), n, n_branches, E, min_sample){
  # evaluate a split at partition branch on covariate j at the splitpoint s
  # index: index of observations in branch
  
  # dividing observation in branch
  index_branch <- index[X[index, j] <= s]
  index_n_branches <- index[X[index, j] > s]
  
  # check wether this split resolves in two reasnable partitions
  if(length(index_branch) < min_sample | length(index_n_branches) < min_sample){
    # remove no longer needed objects from memory
    rm(Q, index, index_branch, index_n_branches)
    return(list('loss' = Inf, j = j, s = s))
  }
  
  # new indicator matrix
  E <- cbind(E, matrix(0, n, 1))
  E[index_n_branches, branch] <- 0
  E[index_n_branches, n_branches] <- 1
  
  # new level estimates
  E_tilde <- Q %*% E
  c_hat <- fastLmPure(E_tilde, Y_tilde)$coefficients

  # resulting new spectral loss
  loss <- loss(Y_tilde, E_tilde %*% c_hat)
  
  # remove no longer needed objects from memory
  rm(Q, Y_tilde, E, E_tilde, c_hat, index, index_branch, index_n_branches)

  return(list('loss' = loss, j = j, s = s))
}

get_all_splitt <- function(branch, X, Y_tilde, Q = diag(length(Y)), n, n_branches, E, multicore = F, min_sample){
  # finds the best splitts for every covariate in branch
  # returns the best splitpoint for every covariate and the resulting loss decrease
  
  # observations belonging to branch
  index <- which(E[, branch] == 1)
  
  # all possible split points
  s <- find_s(X[index, ])
  # itetator for the splits
  iter <- 1:length(s)-1
  
  # evaluate all the relevant splitts
  if(multicore){
      res <- mclapply(iter, function(x)evaluate_splitt(branch, floor(x / dim(s)[1] + 1), 
              s[x + 1], index, X, Y_tilde, Q, n, n_branches, E, min_sample), mc.cores =  n_cores)
  }else{
    res <- lapply(iter, function(x)evaluate_splitt(branch, floor(x / dim(s)[1] + 1), 
                s[x + 1], index, X, Y_tilde, Q, n, n_branches, E, min_sample))
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

predict.SDTree <- function(object, newdata){
  # predict function for the spectral deconfounded tree
  return(predict_outsample(object$tree, newdata))
}

predict.SDforest <- function(object, newdata){
  # predict function for the spectral deconfounded random forest
  # using the mean over all trees as the prediction
  if(is.null(newdata)){
    return(object$f_X_hat)
  }else{
    pred <- do.call(cbind, lapply(object[[2]], function(x){predict_outsample(x$tree, newdata)}))
    return(rowMeans(pred))
  }
}

SDForest <- function(X, Y, nTree, m = 20, cp = 0.005, min_sample = 5, deconfounding = T, mtry = F, multicore = T, pruning = F){
  # Spectral deconfounded random forest using SDTree
  # nTree: number of regression trees for the forest
  # m: maximum number of splits in each tree
  # cp: cost-complexity parameter for each tree
  # min_sample: minimum observations for each partition
  # mtry: number of randomly selected covariates that are available for each split
  # pruning: wether each tree should be automatically pruned using a validation set
  
  # number of observations
  n <- length(Y)
  
  # bootstrap samples
  ind <- lapply(1:nTree, function(x)sample(1:n, n, replace = T))
  
  # estimating all the trees
  if(multicore){
    res <- mclapply(ind, function(i)SDTree(X[i, ], Y[i], m = m, cp = cp, min_sample = min_sample, 
                                           deconfounding = deconfounding, mtry = mtry, fast = T, 
                                           pruning = pruning, Q_type = 1), 
                    mc.cores = n_cores)
  }else{
    res <- lapply(ind, function(i)SDTree(X[i, ], Y[i], m = m, cp = cp, min_sample = min_sample, 
                                         deconfounding = deconfounding, mtry = mtry, fast = T, 
                                         pruning = pruning, Q_type = 1))
  } 
  
  # predict with all trees
  pred <- do.call(cbind, lapply(res, function(x){predict_outsample(x$tree, X)}))
  
  # use mean over trees as final prediction
  f_X_hat <- rowMeans(pred)
  res <- list(f_X_hat = f_X_hat, forest = res)
  class(res) <- 'SDforest'
  return(res)
}