set.seed(24)

n <- 50
x <- rnorm(n)
y <- x**2
a <- factor(sample(letters[1:3], n, replace = TRUE))

expect_error(SDForest(x = x, y = y, envs = a, Q_type = 'no_deconfounding'))
SDForest(x = x, y = y, nTree_leave_out = 2, envs = a, Q_type = 'no_deconfounding')
SDForest(x = x, y = y, nTree_env = 2, envs = a, Q_type = 'no_deconfounding')
ntrees <- c(3, 1, 5)

# need names
expect_error(SDForest(x = x, y = y, nTree_leave_out = ntrees, envs = a, Q_type = 'no_deconfounding'))
expect_error(SDForest(x = x, y = y, nTree_env = ntrees, envs = a, Q_type = 'no_deconfounding'))

names(ntrees) <- levels(a)
fit <- SDForest(x = x, y = y, nTree_leave_out = ntrees, envs = a, Q_type = 'no_deconfounding')
expect_null(fit$nTree_env)
expect_equal(fit$tree_env, ntrees)
expect_equal(length(fit$forest), sum(ntrees))

fit <- SDForest(x = x, y = y, nTree_env = ntrees, envs = a, Q_type = 'no_deconfounding')
expect_null(fit$nTree_leave_out)
expect_equal(fit$tree_env, ntrees)
expect_equal(length(fit$forest), sum(ntrees))