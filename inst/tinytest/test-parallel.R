set.seed(1)

X <- matrix(rnorm(50 * 20), nrow = 50)
Y <- rnorm(50)

SDForest(x = X, y = Y, Q_type = 'no_deconfounding', 
        nTree = 2, mc.cores = 2)
