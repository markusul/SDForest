source("R/SDForest_gpu.r")

library('ranger')

m <- 5
p <- 500
n <- 5000
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

set.seed(7)
data <- simulate_data_nonlinear(1, 100, 100, 1, a = 3)

a <- 3
X_train <- data$X[data$A != a, ]
Y_train <- data$Y[data$A != a]
A_train <- as.matrix(data$A[data$A != a])

# fit SDForest
sdf <- SDForest(x = X_train, y = Y_train, Q_type = 'no_deconfounding')
sdf10 <- SDForest(x = X_train, y = Y_train, A = A_train, gamma = 100, Q_type = 'no_deconfounding')
sdf0 <- SDForest(x = X_train, y = Y_train, A = A_train, gamma = 0, Q_type = 'no_deconfounding')
sdfd <- SDForest(x = X_train, y = Y_train)


sdf2$var_names

A_train

plot(data$X[, data$j], data$Y, pch = data$A, ylim = c(min(data$Y, data$f_X), max(data$Y, data$f_X)))
points(data$X[, data$j], predict(sdf, data.frame(data$X)), col = '#204dca', pch = 20, cex = 0.5)
points(data$X[, data$j], data$f_X, col = '#0c960c', pch = 20, cex = 0.5)


path <- regPath(sdfd, multicore = TRUE, oob = T)
dev.off()
plot(path)
plotOOB(path)

path$varImp_path

plot(path10)
plotOOB(path10)
dev.off()




source("R/SDForest_gpu.r")

library('ranger')

m <- 5
p <- 1000
n <- 1000
q <- 4

data <- simulate_data_nonlinear(q, p, n, m)

X <- data$X
Y <- data$Y


a <- Sys.time()
fit1 <- SDTree(x = X, y = Y)
b <- Sys.time()

source("R/SDForest.r")
c <- Sys.time()
fit2 <- SDTree(x = X, y = Y)
d <- Sys.time()

b - a
d - c

max(fit1$predictions - fit2$predictions)

