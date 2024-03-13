source("R/SDForest.r")
library("GPUmatrix")

m <- 5
p <- 500
n <- 500
q <- 4

data <- simulate_data_nonlinear(q, p, n, m)

X <- data$X

X_gpu <- gpu.matrix(X)

start <- Sys.time()
svd(X)
end <- Sys.time()

print(end - start)

start <- Sys.time()
gpu.svd(X_gpu)
end <- Sys.time()