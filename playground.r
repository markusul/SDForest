source("R/SDForest.r")
library(ggplot2)

p <- 20
n <- 20

set.seed(2024)
data <- simulate_data_nonlinear(20, p, n, 4)


source("R/SDForest.r")
multicore <- T

p <- 500
n <- 500
q <- 20

n_test <- 500

set.seed(2024)
data <- simulate_data_nonlinear(q, p, n + n_test, 4)
data_test <- data
data_test$Y <- data_test$Y[(n+1):(n+n_test)]
data_test$X <- data_test$X[(n+1):(n+n_test),]
data_test$f_X <- data_test$f_X[(n+1):(n+n_test)]

data$X <- data$X[1:n,]
data$Y <- matrix(data$Y[1:n])
data$f_X <- data$f_X[1:n]


  input_data <- data.handler(x = data$X, y = data$Y)
  X <- input_data$X
  Y <- input_data$Y

  # number of observations
  n <- nrow(X)
  # number of covariates
  p <- ncol(X)


  ind <- lapply(1:100, function(x)sample(1:n, n, replace = T))

for(i in ind) svd(X[i, ])$d

m <- cor(X[i, ])

diag(m) <- 0

sort(m, decreasing = T)[1:60]
svd(X[i, ]*1.000000000000001)$d

a <- tryCatch({
  svd(X[i, ])$d
}, error = function(e) {
    warning(paste(e, ':X multipied by number close to 1'))
    return(svd(X[i, ] * 1.000000000000001)$d)})
a
