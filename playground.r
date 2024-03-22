source("R/SDForest.r")

library('ranger')

m <- 5
p <- 500
n <- 500
q <- 4


simulate_data_nonlinear <- function(q, p, n, m, eff = NULL, a = 0){
    #simulate data with confounding and non-linear f_X
    # q: number of confounding covariates in H
    # p: number of covariates in X
    # n: number of observations
    # m: number of covariates with a causal effect on Y

    m <- min(m, p)

    n_A_levels <- a

    if(a == 0){
        A <- rep(0, n)
    }else{
        A_levels <- 1:n_A_levels
        A <- sample(A_levels, n, replace = TRUE)
    }

    alpha_1 <- rnorm(p)
    alpha_2 <- rnorm(q)
    alpha_3 <- rnorm(1)

    # complexity of f_X
    complexity <- 5
    # random parameter for fourier basis
    beta <- runif(m * complexity * 2, -1, 1)

    # random confounding covariates H
    if(q == 0){
        H <- matrix(0, nrow = n, ncol = 1)
        
        # random correlation matrix cov(X, H)
        Gamma <- matrix(0, nrow = 1, ncol = p)

        # random coefficient vector delta
        delta <- 0
    }else{
        H <- matrix(rnorm(n * q, 0, 1), nrow = n)
        H <- H + A %*% t(alpha_2)

        # random correlation matrix cov(X, H)
        Gamma <- matrix(rnorm(q * p, 0, 1), nrow = q)

        # random coefficient vector delta
        delta <- rnorm(q, 0, 1)
    }



    if(!is.null(eff)){
        non_effected <- p - eff
        if(non_effected <= 0) stop('eff must be smaller than p or NULL')
        
        Gamma[, sample(1:p, non_effected)] <- 0
    }



    # random error term
    E <- matrix(rnorm(n * p, 0, 1), nrow = n)

    X <- A %*% t(alpha_1) + H %*% Gamma + E
  
    # random sparse subset of covariates in X
    js <- sample(1:p, m)

    # generate f_X
    f_X <- apply(X, 1, function(x) f_four(x, beta, js))
    
    # generate Y
    Y <- f_X + H %*% delta + A %*% t(alpha_3) + rnorm(n, 0, 0.1)
  
    #return data
    return(list(X = X, Y = Y, f_X = f_X, j = js, beta = beta, H = H, A = A))
}

data <- simulate_data_nonlinear(0, 100, 100, 1, a = 3)

# fit SDForest
sdf <- SDForest(x = data$X, y = data$Y, Q_type = 'no_deconfounding')

a <- unique(data$A)[1]
sdf2 <- SDForest(x = data$X[data$A != a, ], y = data$Y[data$A != a], A = as.matrix(data$A[data$A != a]), gamma = 10, Q_type = 'no_deconfounding')
sdf3 <- SDForest(x = data$X[data$A != a, ], y = data$Y[data$A != a])

sdf2$var_names

plot(data$X[, data$j], data$Y, col = as.factor(data$A), ylim = c(-4, 0))
points(data$X[, data$j], sdf$predictions, col = '#204dca', pch = 20, cex = 0.5)
points(data$X[, data$j], predict(sdf2, data.frame(data$X)), col = '#ae20ca', pch = 20, cex = 0.5)
points(data$X[, data$j], predict(sdf3, data.frame(data$X)), col = '#4a0186', pch = 20, cex = 0.5)
points(data$X[, data$j], data$f_X, col = '#0c960c', pch = 20, cex = 0.5)
