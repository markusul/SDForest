source("R/SDForest.r")

set.seed(42)
data <- simulate_data_nonlinear(q = 20, p = 80, n = 100, m = 2)
X <- data$X
Y <- data$Y
data$j

library(ranger)
fit_ranger <- ranger(Y ~ ., data.frame(X, Y), importance = 'impurity')

# Fit the model
fit <- SDForest(x = X, y = Y, nTree = 100)

imp_ranger <- fit_ranger$variable.importance
imp_sdf <- fit$var_importance

plot(imp_ranger, imp_sdf)


path <- regPath(fit)
plotOOB(path)


plot(path, T)

stablePath <- stabilitySelection(fit)
plot(stablePath, T)

path$cp_min

fit <- prune(fit, cp = path$cp_min)

most_imp <- which.max(fit$var_importance)

dep <- condDependence(fit, most_imp)
plot(dep)


