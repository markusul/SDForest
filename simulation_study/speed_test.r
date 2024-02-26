source("R/SDForest.r")
multicore <- T

p <- 5000
n <- 5000
q <- 20

n_test <- 500

set.seed(2024)
data <- simulate_data_nonlinear(q, p, n + n_test, 4)
data_test <- data
data_test$Y <- data_test$Y[(n+1):(n+n_test)]
data_test$X <- data_test$X[(n+1):(n+n_test),]
data_test$f_X <- data_test$f_X[(n+1):(n+n_test)]

data$X <- data$X[1:n,]
data$Y <- matrix(data$Y[1:n])
data$f_X <- data$f_X[1:n]

start <- Sys.time()
fit <- SDForest(x = data$X, y = data$Y, cp = 0, multicore = multicore)
end <- Sys.time()
print(end - start)

print('fit done')