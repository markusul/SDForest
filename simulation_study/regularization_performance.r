source("R/SDForest.r")
library(parallel)

p <- 500
n <- 500
q <- 20

n_test <- 500
cp_seq <- c(seq(0, 0.1, 0.001), seq(0.1, 0.5, 0.03), seq(0.5, 1, 0.1))

N_rep <- 100

print('start')
res_reg <- mclapply(1:N_rep, function(i){
    data <- simulate_data_nonlinear(q, p, n + n_test, 4)
    data_test <- data
    data_test$Y <- data_test$Y[(n+1):(n+n_test)]
    data_test$X <- data_test$X[(n+1):(n+n_test),]
    data_test$f_X <- data_test$f_X[(n+1):(n+n_test)]
    X_test <- data.frame(data_test$X)
    Y_test <- data_test$Y
    f_X_test <- data_test$f_X

    data$X <- data$X[1:n,]
    data$Y <- matrix(data$Y[1:n])
    data$f_X <- data$f_X[1:n]

    Q <- get_Q(data_test$X, 'trim')

    fit <- SDForest(x = data$X, y = data$Y, cp = 0, multicore = F)

    res <- lapply(cp_seq, function(cp){
        pruned_object <- prune(fit, cp, oob = T)
        pred <- predict(pruned_object, newdata = X_test)
        mse <- mean((pred - Y_test)^2)
        f_mse <- mean((pred - f_X_test)^2)
        sde <- mean((Q%*%(Y_test - pred))^2)
        return(list(mse = mse, f_mse = f_mse, oob.sde = pruned_object$oob_SDloss, 
                    oob.mse = pruned_object$oob_loss, sde = sde, cp = cp))
    })
    res <- do.call(rbind, res)
    res <- matrix(unlist(res), nrow = length(cp_seq), dimnames = list(NULL, colnames(res)))
    return(res)
}, mc.cores = n_cores)
print('done')

save(res_reg, file = "simulation_study/results/regularization_performance.RData")
