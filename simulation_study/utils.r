source("R/SDForest.r")
library(ranger)
library(parallel)

performance_measure <- function(n, p, q, n_test, eff){
    data <- simulate_data_nonlinear(q, p, n+n_test, 4, eff)
    data_train <- data.frame(data$X[1:n,], Y = data$Y[1:n])
    data_test <- data.frame(data$X[(n+1):(n+n_test),], Y = data$Y[(n+1):(n+n_test)])

    fit <- SDForest(Y ~ ., data_train, cp = 0, multicore = T)
    fit2 <- ranger(Y ~ ., data_train, num.trees = 100, importance = 'impurity', mtry = floor(0.9 * ncol(data_train)))

    pred <- predict(fit, data_test)
    pred2 <- predict(fit2, data_test)$predictions

    mse <- (data$f_X[(n+1):(n+n_test)] - pred)
    mse2 <- (data$f_X[(n+1):(n+n_test)] - pred2)

    return(list(SDF = mse, ranger = mse2))
}