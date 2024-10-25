set.seed(1)

n <- 50
X <- matrix(rnorm(n * 20), nrow = n)
Y <- rnorm(n)
A <- rnorm(n)

SDForest(x = X, y = Y, A = A, Q_type = 'no_deconfounding')
SDForest(x = X, y = Y, A = A)
