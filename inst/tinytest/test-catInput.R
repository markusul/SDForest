set.seed(42)

n <- 20
X <- data.frame(y = rnorm(n), 
                x1 = factor(sample(letters[1:3], n, replace = T)), 
                x2 = factor(sample(letters[1:5], n, replace = T), ordered = T), 
                x3 = factor(sample(letters[1:5], n, replace = T), ordered = T), 
                x4 = rnorm(n))

expect_error(SDTree(x = X, y = X$y))
fit <- SDForest(y ~., X, nTree = 2)

fac <- names(X)[sapply(X, function(x){is.factor(x) & !is.ordered(x)})]
expect_equal(ncol(X) + length(levels(X[, fac])) - 3, ncol(fit$X))
expect_true(is.numeric(fit$X))
