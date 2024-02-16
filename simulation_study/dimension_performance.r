source("R/SDForest.r")
library(ranger)
library(parallel)

performance_measure <- function(n, p, q, n_test){
    data <- simulate_data_nonlinear(q, p, n+n_test, 4)
    data_train <- data.frame(data$X[1:n,], Y = data$Y[1:n])
    data_test <- data.frame(data$X[(n+1):(n+n_test),], Y = data$Y[(n+1):(n+n_test)])

    fit <- SDForest(Y ~ ., data_train, cp = 0, multicore = F)
    fit2 <- ranger(Y ~ ., data_train, num.trees = 100, importance = 'impurity', mtry = floor(0.9 * ncol(data_train)))

    pred <- predict(fit, data_test)
    pred2 <- predict(fit2, data_test)$predictions

    mse <- mean((data$f_X[(n+1):(n+n_test)] - pred)^2)
    mse2 <- mean((data$f_X[(n+1):(n+n_test)] - pred2)^2)

    return(c(SDF = mse, ranger = mse2))
}

n <- 500
p <- 500
q <- 20
n_test <- 500

N_rep <- 100

n_seq <- seq(100, 1000, 100)
p_seq <- seq(100, 1000, 100)
q_seq <- seq(0, 100, 10)


print('start')
perf_n <- mclapply(1:N_rep, function(i) sapply(n_seq, function(n) performance_measure(n, p, q, n_test)), mc.cores = n_cores)
save(perf_n, n_seq, file = "simulation_study/results/perf_n.RData")
print('n done')

perf_p <- mclapply(1:N_rep, function(i) sapply(p_seq, function(p) performance_measure(n, p, q, n_test)), mc.cores = n_cores)
save(perf_p, p_seq, file = "simulation_study/results/perf_p.RData")
print('p done')

perf_q <- mclapply(1:N_rep, function(i) sapply(q_seq, function(q) performance_measure(n, p, q, n_test)), mc.cores = n_cores)
save(perf_q, q_seq, file = "simulation_study/results/perf_q.RData")
print('q done')






