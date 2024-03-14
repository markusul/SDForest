source("R/SDForest.r")
library("GPUmatrix")
library('ranger')

m <- 5
p <- 500
n <- 500
q <- 4

data <- simulate_data_nonlinear(q, p, n, m)



rfit <- ranger(x = data.frame(data$X), y = data$Y, num.trees = 100)
fit <- SDForest(x = data$X, y = data$Y, max_size = 100, Q_type = 'no_deconfounding', nTree = 100)

var(data$Y)
rfit$prediction.error
fit$oob_loss


X <- data$X

X_gpu <- gpu.matrix(X)

start <- Sys.time()
svd(X)
end <- Sys.time()

print(end - start)

start <- Sys.time()
gpu.svd(X_gpu)
end <- Sys.time()