library(graphics)
library(ranger)
library(xgboost)

str(airquality)
data <- na.omit(airquality)
str(data)
source("SDForest.r")

Q <- get_Q(data[1:80, -1], type = 1)
Q_2 <- t(Q) %*% Q

sd_objective <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")

    grad <- - Q_2 %*% (labels - preds)
    hess <- Q_2
}

fit <- SDTree(Y = data[1:80, 1], X = data[1:80, -1],min_sample = 2, cp = 0, Q_type = 'pca', )
fit$f_X_hat
print(fit$tree, 'value', 's', 'j')
data <- scale(data)


fit <- SDForest(Y = data[1:80, 1], X = data[1:80, -1], m = 100, nTree = 500, min_sample = 2, cp = 0, mtry = 4, deconfounding = T)
fit$f_X_hat
fit2 <- ranger(Ozone ~ ., data = data[1:80, ], num.trees = 100)

dtrain <- xgb.DMatrix(data = as.matrix(data[1:80, -1]), 
                      label = as.matrix(data[1:80, 1]))
fit3 <- xgb.train(data = dtrain, nrounds = 100)

pred <- predict(fit, data[81:111, -1])
pred2 <- predict(fit2, data[81:111, ])
pred3 <- predict(fit3, xgb.DMatrix(data = as.matrix(data[81:111, -1])))

mean((data[81:111, 1] - pred)**2)
mean((data[81:111, 1] - pred2$predictions)**2)
mean((data[81:111, 1] - pred3)**2)

Q_test <- get_Q(data[81:111, -1], type = 1)
mean((Q_test %*% data[81:111, 1] - Q_test %*% pred)**2)
mean((Q_test %*% data[81:111, 1] - Q_test %*% pred2$predictions)**2)
mean((Q_test %*% data[81:111, 1] - Q_test %*%pred3)**2)


?ranger

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

simulate_data_nonlinear_E <- function(q, p, n, m){
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

    b <- rbinom(n, 1, 0.5)
    e <- rnorm(n, -2, 1)
    e[b == 1] <- rnorm(sum(b), 2, 1)
    E[, 1] <- e

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

data <- simulate_data_nonlinear_E(1, 20, 100, 5)
data

plot(svd(scale(data$X))$d)


plot(density(rnorm(100, 0, 1)))
b <- rbinom(100, 1, 0.5)
e <- rnorm(100, -2, 1)
e[b == 1] <- rnorm(sum(b), 2, 1)
e

plot(density(rnorm(100, -0.5, 1) + rnorm(100, 0.5, 1)))
plot(density(e))



simulate_data_nonlinear_D <- function(q, p, n, m, d = 1){
    #simulate data with confounding and non-linear f_X
    # q: number of confounding covariates in H
    # p: number of covariates in X
    # n: number of observations
    # m: number of covariates with a causal effect on Y

    # random confounding covariates H
    H <- matrix(rnorm(n * q, 0, 1), nrow = n)

    # random correlation matrix cov(X, H)
    Gamma <- matrix(rnorm(q * p, 0, 1), nrow = q)
    non_confounded_covs <- sample(1:p, floor((1-d) * p))
    print(non_confounded_covs)
    Gamma[, non_confounded_covs] <- 0


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



data <- simulate_data_nonlinear_D(1, 20, 100, 5, d = 0.1)
data

quantile(sv$d, 0.5)

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
                                print('aa')
                                sv$u %*% diag(d_pca) %*% t(sv$u)
                            },
                            diag(n)) # no_deconfounding
    return(Q)
}

print(paste(c(1, 2, 3, 4)))


X <- matrix(rnorm(100 * 20), nrow = 100)
get_Q(X, 'pa')

paste(names(modes), collapse = ', ')

sv <- svd(X)

type <- 'pc'
q_hat <- 0
!(type %in% names(modes))


modes <- c('trim' = 1, 'DDL_trim' = 2, 'pca' = 3, 'no_deconfounding' = 4)

Q <- switch(modes[type], sv$u %*% diag(D_tilde) %*% t(sv$u), 
                         diag(n) - sv$u %*% diag(1 - D_tilde) %*% t(sv$u), 
                         {
                            d_pca <- sv$d
                            if(q_hat <= 0) stop("the assumed confounding dimension q_hat must be larger than zero")
                            d_pca[1:q_hat] <- 0
                            print('aa')
                            sv$u %*% diag(d_pca) %*% t(sv$u)
                         },
                         diag(n))
