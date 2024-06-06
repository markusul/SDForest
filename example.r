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


X <- matrix(rnorm(50 * 5), nrow = 50)
Y <- 3 * X[, 1]**2 + rnorm(50)

x <- rnorm(100)
y <- sign(x) * 3 + rnorm(100)
model <- SDTree(x = x, y = y, Q_type = 'no_deconfounding')
model <- prune(model, 1)
plot(model)
pd <- partDependence(model, 1, X = x)
plot(pd)


X <- matrix(rnorm(50 * 20), nrow = 50)
Y <- rnorm(50)
tree <- SDTree(x = X, y = Y)
pruned_tree <- prune(copy(tree), 0.2)
print(tree)
pruned_tree

class(tree)
print.SDTree(tree)
predict(tree, data.frame(X))


set.seed(1)
n <- 50
X <- matrix(rnorm(n * 5), nrow = n)
y <- sign(X[, 1]) * 3 + rnorm(n, 0, 5)
cp <- cvSDTree(x = X, y = y, Q_type = 'no_deconfounding')
cp

library(rpart)
fit <- rpart(y ~ ., data.frame(X, y))
fit$cptable

model <- prune(model, 1)
model



paths <- regPath(model)
data(iris)
tree <- SDForest(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width, 
                 iris, nTree = 10)
tree
varImp(tree)

set.seed(1)
n <- 10
X <- matrix(rnorm(n * 5), nrow = n)
y <- sign(X[, 1]) * 3 + sign(X[, 2]) + rnorm(n)
model <- SDForest(x = X, y = y, Q_type = 'no_deconfounding')
model
paths <- regPath(model)
plotOOB(paths)
plot(paths)
plot(paths, plotly = T)
