set.seed(1)

n <- 50
X <- matrix(rnorm(n * 20), nrow = n)
Y <- X[, 1] * rnorm(n)
A <- sample(1:3, n, prob = c(0.3, 0.3, 0.3), replace = TRUE)
A <- as.factor(A)

SDForest(x = X, y = Y, envs = A, nTree_leave_out = 5)
SDForest(x = X, y = Y, envs = A, nTree_env = 5)
SDForest(x = X, y = Y, envs = A, nTree_leave_out = c('1' = 0, '2' = 3, '3' = 8))
