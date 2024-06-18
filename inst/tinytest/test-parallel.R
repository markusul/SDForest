set.seed(1)

X <- matrix(rnorm(50 * 20), nrow = 50)
Y <- rnorm(50)

fit <- SDForest(x = X, y = Y, Q_type = 'no_deconfounding', 
        nTree = 30, mc.cores = 2)
path <- regPath(fit)
plotOOB(path, sqrt_scale = T)

