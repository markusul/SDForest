source("R/SDForest.r")

library('ranger')

m <- 5
p <- 500
n <- 5000
q <- 4

data <- simulate_data_nonlinear(q, p, n, m)

print(var(data$Y))
rfit <- ranger(x = data.frame(data$X), y = data$Y, num.trees = 100)
print(rfit$prediction.error)
fit <- SDForest(x = data$X, y = data$Y, max_size = 50, Q_type = 'no_deconfounding', nTree = 100, multicore = T)
print(fit$oob_loss)

fit <- SDForest(x = data$X, y = data$Y, max_size = 100, Q_type = 'no_deconfounding', nTree = 100, multicore = T)
print(fit$oob_loss)



fit <- SDForest(x = data$X, y = data$Y, max_size = 700, Q_type = 'no_deconfounding', nTree = 100, multicore = T)
print(fit$oob_loss)
