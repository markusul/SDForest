# Load packages
library(glmnet)
library(parallel)
library(rpart)
library(dplyr)
library(tidyr)

library(rpart)
library(ranger)

simulate_data_step <- function(q, p, n, m, make_tree = F){
    # q: number of confounding covariates in H
    # p: number of covariates in X
    # n: number of observations
    # m: number of partitions

    # minimum number of observations for split
    min_sample <- 2

    # random confounding covariates H
    H <- matrix(rnorm(n * q), nrow = n)

    # random correlation matrix cov(X, H)
    Gamma <- matrix(rnorm(q * p), nrow = q)

    # random coefficient vector delta
    delta <- rnorm(q, 0, 10)

    # relevant covariates
    js <- c()

    # generate X
    if(q == 0){
        # no confounding covariates
        X <- matrix(rnorm(n * p), nrow = n)
    }else{
        # confounding covariates
        X <- H %*% Gamma + matrix(rnorm(n * p), nrow = n)
    }

    # generate tree
    if(make_tree){
        tree <- Node$new(name = '1', value = 0)
    }

    # partitions of observations
    index <- list(1:n)
    for (i in 1:m){
        # get number of observations in each partition
        samples_per_part <- unlist(lapply(index, function(x)length(x)))

        # get potential splits (partitions with enough observations)
        potential_splitts <- which(samples_per_part >= min_sample)

        # sample potential partition to split
        # probability of each partition is proportional to number of observations
        if(length(potential_splitts) == 1){
            branch <- potential_splitts
        }else {
            branch <- sample(potential_splitts, 1,
                            prob = samples_per_part[potential_splitts]/sum(samples_per_part[potential_splitts]))
        }

        # sample covariate to split on
        j <- sample(1:p, 1)
        js <- c(j, js)

        # sample split point
        potential_s <- X[index[[branch]], j]
        s <- rbeta(1, 2, 2) * (max(potential_s) - min(potential_s)) + min(potential_s)

        # split partition
        index <- append(index, list(index[[branch]][X[index[[branch]], j] > s]))
        index[[branch]] <- index[[branch]][X[index[[branch]], j] <= s]

        # add split to tree
        if(make_tree){
            if(tree$height == 1){
                leave <- tree
            }else{
                leaves <- tree$leaves
                leave <- leaves[[which(tree$Get('name', filterFun = isLeaf) == branch)]]
            }
            leave$j <- j
            leave$s <- s
            leave$AddChild(branch, value = 0)
            leave$AddChild(i + 1, value = 0)
        }
     }

    # sample means per partition
    f_X_means <- runif(length(index), -50, 50)

    # generate f_X
    f_X <- rep(0, n)
    for (i in 1:length(index)){
        f_X[index[[i]]] <- f_X_means[i]
    }

    # generate Y
    if(q == 0){
        # no confounding covariates
        Y <- f_X + rnorm(n)
    }else{
        # confounding covariates
        Y <- f_X + H %*% delta + rnorm(n)
    }

    # return data
    if(make_tree){
        # add leave values to tree
        for(l in tree$leaves){
            l$value <- f_X_means[as.numeric(l$name)]
        }
        return(list(X = X, Y = Y, f_X = f_X, j = js, tree = tree, index = index))
    }

    # return data
    return(list(X = X, Y = Y, f_X = f_X, j = js, index = index))
}

simulate_data_nonlinear <- function(q, p, n, m){
    #simulate data with confounding and non-linear f_X
    # q: number of confounding covariates in H
    # p: number of covariates in X
    # n: number of observations
    # m: number of covariates with a causal effect on Y

    # random confounding covariates H
    H <- matrix(rnorm(n * q, 0, 1), nrow = n)

    # random correlation matrix cov(X, H)
    Gamma <- matrix(rnorm(q * p, 0, 1), nrow = q)

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

    # complexity of f_X
    complexity <- 5
    # random parameter for fourier basis
    beta <- runif(m * complexity * 2, -1, 1)
    # generate f_X
    f_X <- apply(X, 1, function(x) f_four(x, beta, js))
    
    # generate Y
    Y <- f_X + H %*% delta + rnorm(n, 0, 0.1)
  
    #return data
    return(list(X = X, Y = Y, f_X = f_X, j = js, beta = beta))
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
        do.call(sum, lapply(1:complexity, function(k) beta[beta_ind[1 + (k-1) *2]] * sin(k * 0.2 * x[j]) + beta[beta_ind[2 + (k-1) *2]] * cos(k * 0.2 * x[j])))
        }))
}

estim_tree_model <- function(X, Y, train_ind, cp = 0){
    # estimate f_X using rpart regression tree
    # X: observed covariates
    # Y: observed response
    # train_ind: indices of training data
    # cp: complexity parameter

    data <- data.frame(X = X, response = Y)

    # cross-validation to select cp
    tree <- rpart(response ~ ., data = data[train_ind, ], control = rpart.control(xval = 20, cp = cp, minsplit = 5, minbucket = 2))
    cptable <- tree$cptable
    cptable <- cptable[cptable[, "nsplit"] < min(n, ncol(X)), ]

    # prune tree
    if(cp == 0){
        tree <- prune(tree, cp = cptable[which.min(cptable[, "xerror"]), "CP"])
    }

    # estimate f_X_hat
    f_X_hat <- predict(tree, newdata = data[-train_ind, ])

    # return estimated f_X_hat
    return(list(f_X_hat = f_X_hat, tree = tree))
}

estim_postprocessed_tree_model <- function(X, Y, train_ind, cp = 0){
    # estimate f_X using regression tree and postprocessing
    # X: observed covariates
    # Y: observed response
    # train_ind: indices of training data
    # cp: complexity parameter

    data <- data.frame(X = X, response = Y)

    # cross-validation to select cp
    tree <- rpart(response ~ ., data = data[train_ind, ], control = rpart.control(xval = 20, cp = cp, minsplit = 5, minbucket = 2))
    cptable <- tree$cptable
    cptable <- cptable[cptable[, "nsplit"] < min(n, ncol(X)), ]

    # prune tree
    if(cp == 0){
        tree <- prune(tree, cp = cptable[which.min(cptable[, "xerror"]), "CP"])
    }

    # estimate f_X_hat
    f_X_hat_train <- predict(tree, newdata = data[train_ind, ])
    f_X_hat <- predict(tree, newdata = data[-train_ind, ])

    # postprocess tree
    rpart_levels <- unique(f_X_hat_train)
    m <- length(rpart_levels)

    # indicator matrix for training data
    E_train <- matrix(0, nrow = length(train_ind), ncol = m)
    for(i in 1:m){
        E_train[f_X_hat_train == rpart_levels[i], i] <- 1
    }

    # calculate Trim transform
    Q <- get_Q(X[train_ind, ], 1)

    # deconfound E and Y
    E_tilde <- Q %*% E_train
    Y_tilde <- Q %*% Y[train_ind]

    # estimate c_hat
    c_hat <- fastLmPure(E_tilde, Y_tilde)$coefficients

    # postprocess f_X_hat
    for(i in 1:m){
        f_X_hat[f_X_hat == rpart_levels[i]] <- c_hat[i]
    }

    # return estimated f_X_hat
    return(list(f_X_hat = f_X_hat, tree = tree))
}

estim_SDTree_model <- function(X, Y, train_ind, multicore = F, cp = 0.0001){
    # estimate f_X using regression tree with deconfounding
    # X: observed covariates
    # Y: observed response
    # train_ind: indices of training data
    # multicore: use multicore
    # cp: complexity parameter

    # estimate tree
    res <- SDTree(X[train_ind, ], Y[train_ind], m = 50, cp = cp, min_sample = 5, deconfounding = T, mtry = F, fast = T, pruning = T, multicore = multicore)

    # predict f_X_hat for test data
    f_X_hat <- predict(res, X[-train_ind, ])

    return(list(f_X_hat = f_X_hat, tree = res$tree))
}

estim_SDTree_model_cv <- function(X, Y, train_ind, multicore = F){
    # estimate f_X using regression tree with deconfounding and cross-validation
    # X: observed covariates
    # Y: observed response
    # train_ind: indices of training data
    # multicore: use multicore

    # cross-validation to select cp
    cp_min <- cv.SDTree(X[train_ind, ], Y[train_ind], multicore = multicore)

    # estimate tree
    res <- SDTree(X[train_ind, ], Y[train_ind], m = 100, cp = cp_min, min_sample = 5, deconfounding = T, mtry = F, fast = T, multicore = multicore)

    # predict f_X_hat for test data
    f_X_hat <- predict(res, X[-train_ind, ])

    return(list(f_X_hat = f_X_hat, tree = res$tree))
}

estim_SDForest <- function(X, Y, train_ind, multicore = F){
    # estimate f_X using random forest with deconfounding
    # X: observed covariates
    # Y: observed response
    # train_ind: indices of training data
    # multicore: use multicore

    # number of covariates
    p <- ncol(X)

    # estimate forest
    res <- SDForest(X[train_ind, ], Y[train_ind], m = 100, cp = 0.0001, mtry = ceiling(p*0.9), nTree = 200, multicore = multicore, pruning = F, deconfounding = T)

    # predict f_X_hat for test data
    f_X_hat <- predict(res, X[-train_ind, ])
    return(list(f_X_hat = f_X_hat, forest = res))
}

estim_ranger <- function(X, Y, train_ind){
    # estimate f_X using ranger
    # X: observed covariates
    # Y: observed response
    # train_ind: indices of training data
    
    # number of covariates
    p <- ncol(X)

    ranger_data <- data.frame(X = X, Y = Y)

    # estimate forest
    res <- ranger(Y ~ ., data = ranger_data[train_ind, ], num.trees = 200, mtry = ceiling(p*0.9), min.node.size = 5)

    # predict f_X_hat for test data
    f_X_hat <- predict(res, ranger_data[-train_ind, ])$predictions
    return(list(f_X_hat = f_X_hat, forest = res))
}

simulate_data_linear <- function(q, p, n, m){
    # q: number of confounding covariates in H
    # p: number of covariates in X
    # n: number of observations
    # m: number of non-zero coefficients in beta

    # random confounding covariates H
    H <- matrix(rnorm(n * q, 0, 1), nrow = n)

    # random correlation matrix cov(X, H)
    Gamma <- matrix(rnorm(q * p, 0, 1), nrow = q)

    # random sparse subset of covariates in X
    js <- sample(1:p, m)

    # random coefficient vector beta
    beta <- sample(c(-1, 1), m, replace = TRUE)

    # random coefficient vector delta
    delta <- rnorm(q, 0, 10)
    # generate data
    if(q == 0){
        # no confounding covariates
        X <- matrix(rnorm(n * p, 0, 1), nrow = n)
        X_used <- X[, js]
        f_X <- X_used %*% beta
        Y <- f_X + rnorm(n)
    }else{
        # confounding covariates
        X <- H %*% Gamma + matrix(rnorm(n * p, 0, 1), nrow = n)
        X_used <- X[, js]
        f_X <- X_used %*% beta
        Y <- f_X + H %*% delta + rnorm(n)
    }

    # true beta
    beta_all <- rep(0, p)
    beta_all[js] <- beta
    return(list(X = X, Y = Y, f_X = f_X, j = js, beta = beta_all))
}

estim_linear_model <- function(X, Y, train_ind){
    # estimate f_X using linear regression
    # X: observed covariates
    # Y: observed response
    # train_ind: indices of training data

    # select best lambda using cross-validation
    cv_model <- cv.glmnet(X[train_ind, ], Y[train_ind], alpha = 1, intercept = FALSE)
    best_lambda <- cv_model$lambda.1se

    # estimate beta using best lambda
    best_model <- glmnet(X[train_ind, ], Y[train_ind], alpha = 1, lambda = best_lambda, intercept = FALSE)
    beta <- best_model$beta

    # estimate f_X_hat
    f_X_hat <- X[-train_ind, ] %*% beta

    # return estimated f_X_hat
    return(list(f_X_hat = f_X_hat, beta = beta))
}

estim_deconfounded_linear_model <- function(X, Y, train_ind){
    # estimate f_X using linear regression with deconfounding
    # X: observed covariates
    # Y: observed response
    # train_ind: indices of training data

    # calculate Trim transform
    Q <- get_Q(X[train_ind, ], 1)

    # deconfound X and Y
    X_tilde <- Q %*% X[train_ind, ]
    Y_tilde <- Q %*% Y[train_ind]

    # select best lambda using cross-validation
    cv_model <- cv.glmnet(X_tilde, Y_tilde, alpha = 1, intercept = FALSE)
    best_lambda <- cv_model$lambda.min

    # estimate beta using best lambda
    best_model <- glmnet(X_tilde, Y_tilde, alpha = 1, 
                            lambda = best_lambda, intercept = FALSE)
    beta <- best_model$beta

    # estimate f_X_hat
    f_X_hat <- X[-train_ind, ] %*% beta

    # return estimated f_X_hat
    return(list(f_X_hat = f_X_hat, beta = beta))
}

estim_deconfounded_linear_model_lambda <- function(X, Y, train_ind){
    # estimate f_X using linear regression with deconfounding
    # if n > p, then this function only uses p randomly selected observations for cross-validation
    # the estimation of beta is done with all observations
    # X: observed covariates
    # Y: observed response
    # train_ind: indices of training data
    
    # number of covariates
    p <- ncol(X)

    # number of observations
    n <- length(train_ind)
    
    # ids of observations used for cross-validation
    n_cv <- train_ind

    # if n > p, then only use p randomly selected observations for cross-validation
    if(length(train_ind) > p){
        n_cv <- sample(train_ind, p)
    }

    # calculate Trim transform
    Q_cv <- get_Q(X[n_cv, ], 1)
    Q <- get_Q(X[train_ind, ], 1)

    # deconfound X and Y
    X_tilde <- Q %*% X[train_ind, ]
    Y_tilde <- Q %*% Y[train_ind]

    X_tilde_cv <- Q_cv %*% X[n_cv, ]
    Y_tilde_cv <- Q_cv %*% Y[n_cv]

    # select best lambda using cross-validation

    cv_model <- cv.glmnet(X_tilde_cv, Y_tilde_cv, alpha = 1, intercept = FALSE)
    best_lambda <- cv_model$lambda.min

    # estimate beta using best lambda
    best_model <- glmnet(X_tilde, Y_tilde, alpha = 1, 
                            lambda = best_lambda, intercept = FALSE)
    beta <- best_model$beta

    # estimate f_X_hat
    f_X_hat <- X[-train_ind, ] %*% beta

    # return estimated f_X_hat
    return(list(f_X_hat = f_X_hat, beta = beta, Q = Q))
}

estim_deconfounded_linear_model_cv <- function(X, Y, train_ind, n_cv = 10){
    # estimate f_X using linear regression with deconfounding
    # cross-validation is used to select lambda with deconfounding for each set separately
    # X: observed covariates
    # Y: observed response
    # train_ind: indices of training data

    # test set
    X_test <- X[-train_ind, ]
    Y_test <- Y[-train_ind]

    # training set
    X <- X[train_ind, ]
    Y <- Y[train_ind]
    
    # number of covariates
    p <- ncol(X)

    # number of observations
    n <- length(Y)

    # number of observations per validation set
    len_test <- floor(n / n_cv)

    # list with all the validation sets
    test_ind <- lapply(1:n_cv, function(x)1:len_test + (x - 1) * len_test)

    # lambda values to test
    lambda_ls <- seq(0, 1, 0.01)

    # estimate performance for each lambda and each validation set
    perf <- lapply(test_ind, function(cv_ind){
        # calculate Trim transform
        Q_cv <- get_Q(X[cv_ind, ], 1)
        Q <- get_Q(X[-cv_ind, ], 1)

        # deconfound X and Y
        X_tilde <- Q %*% X[-cv_ind, ]
        Y_tilde <- Q %*% Y[-cv_ind]

        X_tilde_cv <- Q_cv %*% X[cv_ind, ]
        Y_tilde_cv <- Q_cv %*% Y[cv_ind]

        # estimate lasso for each lambda
        res <- lapply(lambda_ls, function(lambda) glmnet(X_tilde, Y_tilde, alpha = 1, intercept = FALSE, lambda = lambda))
        
        # calculate performance for each lambda
        perf <- lapply(res, function(model) sum((predict(model, X_tilde_cv) - Y_tilde_cv)**2) / len_test)
        return(perf)
    })

    perf <- matrix(unlist(perf), ncol = n_cv, byrow = FALSE)
    # mean performance for each lambda over all validation sets
    perf_mean <- apply(perf, 1, mean)

    # lambda with minimum loss
    min_lambda <- which.min(perf_mean)
    lambda.min <- lambda_ls[min_lambda]
    # lambda with minimum loss + 1 standard deviation
    lambda.se <- max(lambda_ls[lambda_ls < lambda.min + sd(perf[min_lambda, ])])

    # calculate Trim transform
    Q <- get_Q(X, 1)
    X_tilde <- Q %*% X
    Y_tilde <- Q %*% Y

    # estimate beta using best lambda
    best_model <- glmnet(X_tilde, Y_tilde, alpha = 1, 
                            lambda = lambda.se, intercept = FALSE)
    beta <- best_model$beta

    # estimate f_X_hat
    f_X_hat <- X_test %*% beta

    # return estimated f_X_hat
    return(list(f_X_hat = f_X_hat, beta = beta, Q = Q))
}

performance <- function(f_X, f_X_hat){
    # measure squared error of f_X_hat
    # f_X: true f_X
    # f_X_hat: estimated f_X_hat

    # calculate squared error
    se <- sum((f_X - f_X_hat)^2) / length(f_X)
    return(se)
}

performance_comparison <- function(q, p, n, m, data_gen, estimators){
    # compare performance of estimators given data generation function
    # q: number of confounding covariates in H
    # p: number of covariates in X
    # n: number of observations
    # m: number of non-zero coefficients in beta
    # data_gen: data generation function
    # estimators: list of estimators

    # generate data
    data <- data_gen(q, p, n, m)

    # training data
    train_ind <- 101:n

    # estimate f_X_hat for each estimator
    # and calculate performance
    perf <- list()
    for(estimator in estimators){
        f_X_hat <- estimator(data$X, data$Y, train_ind)$f_X_hat
        perf <- append(perf, performance(data$f_X[-train_ind], f_X_hat))
    }
    names(perf) <- names(estimators)

    # return performance
    return(perf)
}

aggregated_performance_comparison <- function(N, q, p, n, m, data_gen, estimators){
    # compare aggregated performance of estimators given data generation function
    # N: number of simulations
    # q: number of confounding covariates in H
    # p: number of covariates in X
    # n: number of observations
    # m: number of non-zero coefficients in beta
    # data_gen: data generation function
    # estimators: list of estimators

    # run performance comparison for each simulation
    res <- mclapply(1:N, function(x)performance_comparison(q, p, n, m, 
                    data_gen, estimators), mc.cores = n_cores)
    res <- unlist(res)
    res <- data.frame('error' = res, 'estimator' = names(res))

    # aggregate performance over simulations
    res <- res %>% group_by(estimator) %>% summarise(mean = mean(error), 
                                                    q_05 = quantile(error, 0.05), 
                                                    q_95 = quantile(error, 0.95))
    return(res)
}

N_performance_comparison <- function(N, q, p, n, m, data_gen, estimators){
    # N: number of simulations
    # q: number of confounding covariates in H
    # p: number of covariates in X
    # n: number of observations
    # m: number of non-zero coefficients in beta
    # data_gen: data generation function
    # estimators: list of estimators

    # run performance comparison for each simulation
    res <- mclapply(1:N, function(x)performance_comparison(q, p, n, m, 
                    data_gen, estimators), mc.cores = n_cores)
    res <- unlist(res)
    res <- data.frame('error' = res, 'estimator' = names(res), 'q' = q, 'p' = p, 'n' = n)

    # return aggregated performance
    return(res)
}