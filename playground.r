source("R/SDForest.r")

library('ranger')

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
    dep <- list(alpha_1 = alpha_1, alpha_2 = alpha_2, alpha_3 = alpha_3, beta = beta, delta = delta, Gamma = Gamma)
    return(list(X = X, Y = Y, f_X = f_X, j = js, beta = beta, H = H, A = A, dep = dep))
}

set.seed(2)
data <- simulate_data_nonlinear(1, 10, 300, 1, a = 3)

n <- 2000
q <- 1
p <- 10

a_seq <- seq(-10, 10, 1)
a_seq <- 3
mse_list <- c()
mse_inv_list <- c()

plot(data$X[, data$j], data$Y, ylim = c(-10, 10), xlim = c(-12, 12))
points(data$X[, data$j], data$f_X, col = 'green', pch = 20, cex = 0.5)
for(i in a_seq){
    A <- rep(i, n)
    H <- matrix(rnorm(n * q, 0, 1), nrow = n)
    H <- H + A %*% t(data$dep$alpha_2)
    E <- matrix(rnorm(n * p, 0, 1), nrow = n)
    X <- A %*% t(data$dep$alpha_1) + H %*% data$dep$Gamma + E
    f_X <- apply(X, 1, function(x) f_four(x, data$beta, data$j))

    y_a <- f_X + H %*% data$dep$delta + A %*% t(data$dep$alpha_3) + rnorm(n_a, 0, 0.1)

    y_pred_inv <- predict(sdf10, data.frame(X))
    y_pred <- predict(sdf, data.frame(X))
    mse_inv <- quantile(abs(y_pred_inv - y_a), 0.99)
    mse <- quantile(abs(y_pred - y_a), 0.99)

    mse_list <- c(mse_list, mse)
    mse_inv_list <- c(mse_inv_list, mse_inv)

    points(X[, data$j], y_a, col = i + 4, pch = 20, cex = 0.5)
}

plot(a_seq, mse_list, col = 'black')
points(a_seq, mse_inv_list, col = 'red', pch = 20, cex = 0.5)
legend('topleft', legend = c('gam = 1', 'gam = 100'), col = c('black', 'red'), pch = c(1, 20), cex = 1)

points(data$X[, data$j], predict(sdf, data.frame(data$X)), col = '#204dca', pch = 20, cex = 0.5)
points(data$X[, data$j], predict(sdf10, data.frame(data$X)), col = '#aa0caf', pch = 20, cex = 0.5)
points(data$X[, data$j], predict(sdf0, data.frame(data$X)), col = '#af0c58', pch = 20, cex = 0.5)
points(data$X[, data$j], predict(sdfd, data.frame(data$X)), col = '#0c9caf', pch = 20, cex = 0.5)
points(data$X[, data$j], predict(sdfda, data.frame(data$X)), col = '#84e64c', pch = 20, cex = 0.5)





a <- 3
X_train <- data$X[data$A != a, ]
Y_train <- data$Y[data$A != a]
A_train <- as.matrix(data$A[data$A != a])

# fit SDForest
sdf <- SDForest(x = X_train, y = Y_train, Q_type = 'no_deconfounding')
sdf10 <- SDForest(x = X_train, y = Y_train, A = A_train, gamma = 10, Q_type = 'no_deconfounding')
sdf0 <- SDForest(x = X_train, y = Y_train, A = A_train, gamma = 0, 
    Q_type = 'no_deconfounding')
sdfd <- SDForest(x = X_train, y = Y_train)
sdfda <- SDForest(x = X_train, y = Y_train, A = A_train, gamma = 0)



plot(data$X[, data$j], data$Y, pch = data$A, ylim = c(min(data$Y, data$f_X), max(data$Y, data$f_X)))
points(data$X[, data$j], predict(sdf, data.frame(data$X)), col = '#204dca', pch = 20, cex = 0.5)
points(data$X[, data$j], predict(sdf10, data.frame(data$X)), col = '#aa0caf', pch = 20, cex = 0.5)
points(data$X[, data$j], predict(sdf0, data.frame(data$X)), col = '#af0c58', pch = 20, cex = 0.5)
points(data$X[, data$j], predict(sdfd, data.frame(data$X)), col = '#0c9caf', pch = 20, cex = 0.5)
points(data$X[, data$j], predict(sdfda, data.frame(data$X)), col = '#84e64c', pch = 20, cex = 0.5)
points(data$X[, data$j], data$f_X, col = '#0c960c', pch = 20, cex = 0.5)
title('n = 200, p = 200, q = 1')
legend('bottomright', legend = c('gam = 1', 'gam = 100', 'gam = 0', 'trim', 'gam = 0 + trim', 'true causal'), 
    col = c('#204dca', '#aa0caf', '#af0c58', '#0c9caf', '#84e64c', '#0c960c'), pch = 20, cex = 1)

path <- regPath(sdfd, multicore = TRUE, oob = T)
dev.off()
plot(path)
plotOOB(path)

path$varImp_path

plot(path10)
plotOOB(path10)
dev.off()


n_train <- 300
r <- 2
q <- 1
p <- 10

A <- matrix(rnorm(n_train * r, 0, 1), nrow = n_train)
H <- matrix(rnorm(n_train * q, 0, 0.01), nrow = n_train)

gamma <- matrix(rnorm(r * p, 0, 1), nrow = r)
delta <- matrix(rnorm(q * p, 0, 1), nrow = q)
X <- A %*% gamma + H %*% delta + matrix(rnorm(n_train * p, 0, 1), nrow = n_train)
Y <- 3 * X[, 2] + 3 * X[, 3] + H - 0.1 * A[, 1] + rnorm(n_train, 0, 1)

n_test <- 2000
A_test <- matrix(rnorm(n_test * r, 0, 1), nrow = n_test) * sqrt(10)
H_test <- matrix(rnorm(n_test * q, 0, 0.1), nrow = n_test)
X_test <- A_test %*% gamma + H_test %*% delta + matrix(rnorm(n_test * p, 0, 1), nrow = n_test)
Y_test <- 3 * X_test[, 2] + 3 * X_test[, 3] + H_test - 0.1 * A[, 1] + rnorm(n_test, 0, 1)


W <- get_W(A, gamma = 10)
fit_linA <- lm.fit(x = W %*% X, y = W %*% Y)$coefficients
fit_lin <- lm.fit(x = X, y = Y)$coefficients

pred_lin <- X_test %*% fit_lin
perf_lin <- quantile(abs(Y_test - pred_lin), seq(0, 1, 0.05))

pred_linA <- X_test %*% fit_linA
perf_linA <- quantile(abs(Y_test - pred_linA), seq(0, 1, 0.05))

plot(perf_lin, ylim = c(0, max(perf_lin, perf_linA)))
points(perf_linA, col = 'red')

plot(X_test[, 2], Y_test)
points(X_test[, 2], pred_lin, col = 'blue', pch = 20, cex = 0.5)
points(X_test[, 2], pred_linA, col = 'red', pch = 20, cex = 0.5)
points(X_test[, 2], 3 * X_test[, 2] + 3 * X_test[, 3], col = 'green', pch = 20, cex = 0.5)

n <- 300
H <- rnorm(n)
A <- (rbinom(n, 1, 0.5) - 1) * 2
X <- A + H + rnorm(n, 0, 1)
Y <- X + 2 * H + rnorm(n, 0, 1)

H_test <- rnorm(n)
A_test <- -4
X_test <- A_test + H_test + rnorm(n, 0, 1)
Y_test <- X_test + 2 * H_test + rnorm(n, 0, 1)


W <- get_W(matrix(A), gamma = 1000)
fit_linA <- lm.fit(x = W %*% X, y = W %*% Y)$coefficients
fit_lin <- lm.fit(x = matrix(X), y = Y)$coefficients

pred_lin <- X_test * fit_lin
pred_linA <- X_test * fit_linA


plot(X, Y, ylim = c(-20, 20), xlim = c(-10, 10))
points(X_test, Y_test, col = 'red')
points(X_test, pred_lin, col = 'blue', pch = 20, cex = 0.5)
points(X_test, pred_linA, col = 'green', pch = 20, cex = 0.5)

source("R/SDForest_gpu.r")

library('ranger')

m <- 5
p <- 5000
n <- 5000
q <- 4

data <- simulate_data_nonlinear(q, p, n, m)

X <- data$X
Y <- data$Y

a <- Sys.time()
fit1 <- SDTree(x = X, y = Y, cp = 0.01, Q_type = 'no_deconfounding')
b <- Sys.time()
b - a

#source("R/SDForest.r")
c <- Sys.time()
fit2 <- rpart::rpart(Y ~ ., data = data.frame(X, Y), cp = 0.01)
d <- Sys.time()

b - a
d - c
fit2$y == Y
predict(fit2, data.frame(X)) - fit1$predictions
data$j
sort(fit1$var_importance, decreasing = TRUE)[1:6]


data <- simulate_data_nonlinear(q, 3000, 3000, m)
X <- data$X
Y <- data$Y

X_gpu <- gpu.matrix(X)
Y_gpu <- gpu.matrix(Y)


c <- Sys.time()
lj <- lapply(1:10000, function(i)t(X_gpu) %*% X_gpu)
d <- Sys.time()

b - a
d - c



library(GPUmatrix)



tot <- 5000 * 5000
 
2e+07
1000000000
800000000
k <- tot / 1000

X <- matrix(0, nrow = k, ncol = 1000)
X_gpu <- gpu.matrix(X)

X_gpu %*% t(X_gpu)

rm(X)
gc()
X_gpu <- as.matrix(X_gpu)

rm(X_gpu)


for(i in 1:100){
    print(i)
    X_gpu <- gpu.matrix(X)
    Y <- t(X_gpu) %*% X
}

torch::