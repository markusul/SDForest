library(SDForest)

set.seed(42)
sim_data <- simulate_data_nonlinear(q = 2, p = 150, n = 100, m = 2)
X <- sim_data$X
Y <- sim_data$Y
train_data <- data.frame(X, Y)
sim_data$j


tree_plain_cv <- cvSDTree(Y ~ ., train_data, Q_type = "no_deconfounding")
tree_plain <- SDTree(Y ~ ., train_data, Q_type = "no_deconfounding", cp = 0)

tree_causal_cv <- cvSDTree(Y ~ ., train_data)
tree_causal <- SDTree(Y ~ ., train_data, cp = 0)

path_plain <- regPath(tree_plain)
plot(path_plain, T)

path <- regPath(tree_causal)
plot(path, T)

tree_plain <- prune(tree_plain, cp = tree_plain_cv$cp_min)
tree_causal <- prune(tree_causal, cp = tree_causal_cv$cp_min)
plot(tree_causal)
plot(tree_plain)



fit_ranger <- ranger::ranger(Y ~ ., train_data, importance = 'impurity')

fit <- SDForest(x = X, y = Y, nTree = 10, Q_type = 'pca', q_hat = 2)
fit <- SDForest(Y ~ ., train_data)
fit

imp_ranger <- fit_ranger$variable.importance
imp_sdf <- fit$var_importance
imp_col <- rep('black', length(imp_ranger))
imp_col[sim_data$j] <- 'red'

plot(imp_ranger, imp_sdf, col = imp_col, pch = 20,
     xlab = 'ranger', ylab = 'SDForest', 
     main = 'Variable Importance')

path <- regPath(fit)
plotOOB(path)
plot(path, plotly = T)

stablePath <- stabilitySelection(fit)
plot(stablePath, plotly = T)

fit <- prune(fit, cp = path$cp_min)

dep <- partDependence(fit, data$j[1])
plot(dep, n_examples = 100)







x <- rnorm(100)
y <- sign(x) * 3 + rnorm(100)
model <- SDTree(x = x, y = y, Q_type = 'no_deconfounding')
plot(model)
model

model$predictions
print(model)
predict(model, newdata = data.frame(X = x))

library(tinytest)
library(SDForest)




library(SDForest)

set.seed(42)
# simulation of confounded data
sim_data <- simulate_data_nonlinear(q = 2, p = 150, n = 100, m = 2)
X <- sim_data$X
Y <- sim_data$Y
train_data <- data.frame(X, Y)
# causal parents
sim_data$j

fit <- SDForest(Y ~ ., train_data)
fit


causal_Tree <- SDTree(Y ~ ., train_data, cp = 0)

# plot the causal tree
plot(causal_Tree)

