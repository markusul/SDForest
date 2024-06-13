set.seed(1)

X <- matrix(rnorm(50 * 20), nrow = 50)
Y <- rnorm(50)

SDForest(x = X, y = Y, Q_type = 'no_deconfounding', 
         mc.cores = 2, nTree = 2)