source("R/SDForest.r")
library(rpart)
library(rattle)

set.seed(42)
data <- simulate_data_nonlinear(q = 2, p = 150, n = 100, m = 2)
X <- data$X
Y <- data$Y
data$j


tree_rpart <- rpart(Y ~ ., data.frame(X, Y), method = 'anova', cp = 0, minsplit = 1, minbucket = 5)
cp_min <- tree_rpart$cptable[which.min(tree_rpart$cptable[, 'xerror']), 'CP']

tree_rpart <- rpart::prune(tree_rpart, cp = cp_min)
fancyRpartPlot(tree_rpart)

tree_plain <- SDTree(Y ~ ., data.frame(X, Y), Q_type = "no_deconfounding", cp = cp_min)
plot(tree_plain)

tree_cv <- cv.SDTree(Y ~ ., data.frame(X, Y))
tree_causal <- SDTree(Y ~ ., data.frame(X, Y), cp = 0)

path <- regPath(tree_causal)
plot(path)

tree_causal <- prune(tree_causal, cp = tree_cv$cp_min)
plot(tree_causal)

library(ranger)
fit_ranger <- ranger(Y ~ ., data.frame(X, Y), importance = 'impurity')

# Fit the model
fit <- SDForest(x = X, y = Y)
fit

imp_ranger <- fit_ranger$variable.importance
imp_sdf <- fit$var_importance

plot(imp_ranger, imp_sdf)


path <- regPath(fit)
plotOOB(path)

path

plot(path, plotly = T)

stablePath <- stabilitySelection(fit)
plot(stablePath, plotly = T)

path$cp_min

fit <- prune(fit, cp = path$cp_min)

most_imp <- which.max(fit$var_importance)

dep <- partDependence(fit, most_imp)
plot(dep)



dep <- partDependence(tree_causal, most_imp, X = data$X)
plot(dep)






x <- rnorm(100)
y <- sign(x) * 3 + rnorm(100)
model <- SDTree(x = x, y = y, Q_type = 'no_deconfounding')
plot(model)
model

model$predictions
print(model)
predict(model, newdata = data.frame(X = x))


library(SDForest)

set.seed(1)

X <- matrix(rnorm(50 * 20), nrow = 50)
Y <- rnorm(50)

tree <- SDTree(x = X, y = Y, Q_type = 'no_deconfounding', 
               cp = 0, min_sample = 5)
min(table(tree$predictions)) >= 5
rpart_tree <- rpart::rpart(y ~ ., data.frame(y = Y, X), 
                    control = rpart::rpart.control(minbucket = 5, 
                                            cp = 0, 
                                            minsplit = 10, 
                                            xval = 0))
pruned_tree <- prune(copy(tree), 0.1)
pruned_rpart_tree <- rpart::prune(rpart_tree, 0.1)

expect_equal(tree$predictions, predict(tree, data.frame(X)))
expect_false(all(tree$predictions == predict(pruned_tree, data.frame(X))))


expect_equal(tree$predictions, as.vector(predict(rpart_tree)))
expect_equal(predict(pruned_tree, data.frame(X)), as.vector(predict(pruned_rpart_tree)))
