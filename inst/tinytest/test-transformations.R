set.seed(1)
#get_Q

# no deconfounding
X <- matrix(rnorm(50 * 20), nrow = 50)
Q_plain <- get_Q(X, 'no_deconfounding')
expect_equal(Q_plain %*% X, X)

#dimensions
n <- 50
X <- matrix(rnorm(n * 20), nrow = n)
Q_trim <- get_Q(X, "trim")
Q_pca <- get_Q(X, "pca", q_hat = 5)
Q_plain <- get_Q(X, "no_deconfounding")

# Q_trim
expect_equal(n, nrow(Q_trim))
expect_equal(n, ncol(Q_trim))

# Q_pca
expect_equal(n, nrow(Q_pca))
expect_equal(n, ncol(Q_pca))

# Q_plain
expect_equal(n, nrow(Q_plain))
expect_equal(n, ncol(Q_plain))


#dimension 1
n <- 50
X <- matrix(rnorm(n * 1), nrow = n)
expect_warning(get_Q(X, "trim"))

#atent > 0
n <- 50
X <- matrix(rnorm(n * 2), nrow = n)
expect_error(get_Q(X, "pca"))

#get_W
#no w
X <- matrix(rnorm(50 * 20), nrow = 50)
W_plain <- get_W(X, gamma = 1)
expect_equal(W_plain %*% X, X)

#dimensions
n <- 50
X <- matrix(rnorm(n * 20), nrow = n)
W <- get_W(X, gamma = 0.5)
# W
expect_equal(n, nrow(W))
expect_equal(n, ncol(W))

#regression residuals
n <- 50
X <- matrix(rnorm(n), nrow = n)
Y <- 3 * X + rnorm(n)
W <- get_W(X, gamma = 0)
expect_equal(lm.fit(y = Y, x = X)$residuals, as.vector(W %*% Y))
