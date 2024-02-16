simulate_data_nonlinear <- function(q, p, n, m){
    #simulate data with confounding and non-linear f_X
    # q: number of confounding covariates in H
    # p: number of covariates in X
    # n: number of observations
    # m: number of covariates with a causal effect on Y

    # random confounding covariates H
    H <- matrix(rnorm(n * q, 0, 1), nrow = n)

    # random correlation matrix cov(X, H)
    Gamma <- matrix(rnorm(q * p, 0, 1), nrow = q)

    # random coefficient vector delta
    delta <- rnorm(q, 0, 1)

    # random error term
    E <- matrix(rnorm(n * p, 0, 1), nrow = n)

    if(q == 0){
        X <- E
    }else{
        X <- H %*% Gamma + E
    }
  
    # random sparse subset of covariates in X
    js <- sample(1:p, m)

    # complexity of f_X
    complexity <- 5
    # random parameter for fourier basis
    beta <- runif(m * complexity * 2, -1, 1)
    # generate f_X
    f_X <- apply(X, 1, function(x) f_four(x, beta, js))
    
    # generate Y
    Y <- f_X + H %*% delta + rnorm(n, 0, 0.1)
  
    #return data
    return(list(X = X, Y = Y, f_X = f_X, j = js, beta = beta))
}

f_four <- function(x, beta, js){
    # function to generate f_X
    # x: covariates
    # beta: parameter vector
    # js: relevant covariates

    # number of relevant covariates
    m <- length(js)

    # complexity of f_X
    complexity <- length(beta) / (2 * m)

    # calculate f_X
    do.call(sum, lapply(1:m, function(i) {
        j <- js[i]
        # select beta for covariate j
        beta_ind <- 1:(2*complexity) + (i-1) * 2 * complexity
        # calculate f_X_j
        do.call(sum, lapply(1:complexity, function(k) beta[beta_ind[1 + (k-1) *2]] * sin(k * 0.1 * x[j]) + beta[beta_ind[2 + (k-1) *2]] * cos(k * 0.1 * x[j])))
        }))
}

true_function <- function(beta, js){
    res <- list(beta = beta, js = js)
    class(res) <- 'true_function'
    return(res)
}

predict.true_function <- function(object, newdata){
    f_X <- apply(newdata, 1, function(x) f_four(x, object$beta, object$js))
    return(f_X)
}

plotDep <- function(object, n_examples = 19){
  ggdep <- ggplot2::ggplot() + ggplot2::theme_bw()
  preds <- object$preds
  x_seq <- object$x_seq
  
  sample_examples <- sample(1:dim(preds)[2], n_examples)
  for(i in sample_examples){
      pred_data <- data.frame(x = x_seq, y = preds[, i])
      ggdep <- ggdep + ggplot2::geom_line(data = pred_data, ggplot2::aes(x = x, y = y), col = 'grey')
  }

  ggdep <- ggdep + ggplot2::geom_line(data = data.frame(x = x_seq, y = object$preds_mean), 
                    ggplot2::aes(x = x, y = y), col = '#08cbba', linewidth = 1.5)
  ggdep <- ggdep + ggplot2::geom_point(data = data.frame(x = object$xj, y = -5), 
                    ggplot2::aes(x = x, y = y), col = 'black', size = 1,shape = 108)
  ggdep <- ggdep + ggplot2::ylab('f(x)') + ggplot2::ggtitle('Conditional dependence')
  ggdep <- ggdep + ggplot2::xlim(quantile(object$xj, 0.05), quantile(object$xj, 0.95))
  if(is.character(object$j)){
    ggdep <- ggdep + ggplot2::xlab(object$j)
  }else{
    ggdep <- ggdep + ggplot2::xlab(paste('x', object$j, sep = ''))
  }
  ggdep + ggplot2::ylim(-5, 5)
}

source("R/SDForest.r")
library(ggplot2)

p <- 40
n <- 40

set.seed(2024)
data <- simulate_data_nonlinear(20, p, n, 4)

true_f <- true_function(data$beta, data$j)

dep_f_1 <- condDependence(true_f, data$j[1], data$X)
dep_f_2 <- condDependence(true_f, data$j[2], data$X)
dep_f_3 <- condDependence(true_f, data$j[3], data$X)
dep_f_4 <- condDependence(true_f, data$j[4], data$X)

set.seed(2024)
gridExtra::grid.arrange(plotDep(dep_f_1), plotDep(dep_f_2), plotDep(dep_f_3), plotDep(dep_f_4), ncol = 2)


dat <- data.frame(X = data$X, Y = data$Y)
dat <- data.frame(dat)


fit <- SDForest(x = data$X, y = data$Y, cp = 0.1, nTree = 3)
fit$var_importance


boxplot(colMeans(dat))

fit <- SDForest(Y ~ ., data = dat, cp = 0, multicore = T)
data$j

most_important <- sort(varImp(fit), decreasing = T)[1]
sort(varImp(fit), decreasing = T)[1:6]
names(most_important)

plot(fit$var_importance)

reg_path <- regPath(fit, oob = T)
plot(reg_path, T)

stable_path <- stabilitySelection(fit)
plot(stable_path, T)

reg_path$loss_path



plotOOB(reg_path)


dep_1 <- condDependence(fit, dat, data$j[1])
dep_2 <- condDependence(fit, dat, data$j[2])
dep_3 <- condDependence(fit, dat, data$j[3])
dep_4 <- condDependence(fit, dat, data$j[4])

gridExtra::grid.arrange(plotDep(dep_1), plotDep(dep_2), plotDep(dep_3), plotDep(dep_4), ncol = 2)



library(ranger)
fit2 <- ranger(Y ~ ., data = dat, num.trees = 100, importance = 'impurity', mtry = 90)
plot(fit2$variable.importance)

ranger_fun <- function(object){
    res <- list(model = object)
    class(res) <- 'ranger_fun'
    return(res)
}
predict.ranger_fun <- function(object, newdata){
    return(predict(object$model, newdata)$predictions)
}

ranger_fit <- ranger_fun(fit2)

dep_r_1 <- condDependence(ranger_fit, dat, data$j[1])
dep_r_2 <- condDependence(ranger_fit, dat, data$j[2])
dep_r_3 <- condDependence(ranger_fit, dat, data$j[3])
dep_r_4 <- condDependence(ranger_fit, dat, data$j[4])

gridExtra::grid.arrange(plotDep(dep_r_1), plotDep(dep_r_2), plotDep(dep_r_3), plotDep(dep_r_4), ncol = 2)


data <- simulate_data_nonlinear(0, 1, 100, 1)
dat <- data.frame(X = data$X, Y = data$Y)
dat <- data.frame(dat)

fit <- SDTree(Y ~ ., data = dat, cp = 0, Q_type = 'no_deconfounding')


plot(data$X, data$Y)
points(data$X, fit$predictions, col = 'red', pch = 20)

path <- regPath(fit)
plot(path)


res <- SDTree(Ozone ~ ., data = airquality, cp = 0)
res$var_importance
res$tree$Get('cp')
prune(res, 1)$var_importance

start_time <- Sys.time()
a <- regPath(res)
print(Sys.time() - start_time)

rm(a)


a$loss_path

prune(res, 0.1)$var_importance

res$oob_SDloss

start <- Sys.time()
res <- SDForest(Ozone ~ ., data = airquality, nTree = 100, cp = 0)
print(Sys.time() - start)


cp_seq <- c(seq(0, 0.1, 0.001), seq(0.1, 0.5, 0.02), seq(0.5, 1, 0.05))
length(cp_seq)


class(a)
print(a)
a$cp_min
plot(a)

res$forest
dep <- condDependence(res, na.omit(airquality), 'Temp')
plot(dep)

predict(res, na.omit(airquality))

library(dplyr)
library(tidyr)
imp_data <- data.frame(a$varImp_path, cp = a$cp)
imp_data <- gather(imp_data, key = 'covariate', value = 'importance', -cp)
imp_data

library(ggplot2)
ggplot(imp_data, aes(x = cp, y = importance, col = covariate)) +
    geom_line() + 
    theme_bw()


loss_data <- data.frame(a$loss_path, cp = a$cp)
gg_sde <- ggplot(loss_data, aes(x = cp, y = oob.SDE)) +
    geom_line() + 
    theme_bw()

gg_mse <- ggplot(loss_data, aes(x = cp, y = oob.MSE)) +
    geom_line() + 
    theme_bw()

gg_sde
gg_mse

res$forest[[1]]$var_importance
res$var_importance
a <- prune(res, 1)

a$forest[[1]]$var_importance
a$var_importance

res$oob_SDloss
a$oob_SDloss


res$forest
a <- lapply(res$forest, function(tree){prune_tree(tree, 0.1)})
a


prune_tree(res, 0.01)$predictions


min(res$tree$Get('cp'), na.rm = T)

res$var_importance == varImp(res$tree, 5)

varImp(res)

res$var_names

res <- SDForest(Ozone ~ ., data = airquality)
res$var_importance == varImp(res)


plot(b$predictions)
points(dat$Y, col = 'red', pch = 20)
points(data$f_X, col = 'green', pch = 20)
b

f <- condDependence(b, dat, data$j[1])
plot(f)

plot(a$var_importance)



library(ranger)
c <- ranger(Y ~ ., data = dat, num.trees = 100, importance = 'impurity')

dep_ranger <- condDependence(c, dat, data$j[1])


ggdep <- plot(f)
ggdep + geom_line(data = data.frame(x = dep_f$x_seq, y = dep_f$preds_mean), aes(x = x, y = y), col = 'red')



plot(data$X[, data$j[1]], data$Y)
points(data$X[, data$j[1]], a$predictions, col = 'red', pch = 20)
points(data$X[, data$j[1]], b$predictions, col = 'blue', pch = 20)
points(data$X[, data$j[1]], data$f_X, col = 'green', pch = 20)

library(ranger)
c <- ranger(Y ~ ., data = dat, num.trees = 100, importance = 'impurity')
points(data$X[, data$j[1]], c$predictions, col = 'purple', pch = 2, cex = 0.5)



plot(c$variable.importance/max(c$variable.importance), ylim = c(0, 1))
points(b$var_imp/max(b$var_imp), col = 'red', pch = 3)
points(a$var_imp/max(a$var_imp), col = 'blue', pch = 4)
points(rep(1, length(data$j)), x = data$j, col = 'green', pch = 2)
grid()


sort(b$var_importance, decreasing = T)[1:6]
sort(a$var_importance, decreasing = T)[1:6]
sort(c$variable.importance, decreasing = T)[1:6]
data$j

start_time <- Sys.time()
suppressWarnings({
a <- lapply(1:40, function(i) SDTree(Y ~ ., dat[sample(1:n, n, replace = T), ], Q_type = 'DDL_trim', cp = 0, max_leaves = 400))
})
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
a <- SDForest(Y ~ ., dat, Q_type = 'DDL_trim', multicore = T, mtry = p, cp = 0, max_leaves = 400, nTree = 40)
end_time <- Sys.time()
end_time - start_time


a$forest[[3]]$tree$height

tree_heights <- unlist(lapply(a$forest, function(x)x$tree$height))
mean(tree_heights)
min(tree_heights)
max(tree_heights)

nTree <- 10
n <- nrow(dat)
max_leaves <- 300
cp <- 0
min_sample <- 5
Q <- get_Q(dat[, - ncol(dat)], 'DDL_trim')
mtry <- floor(sqrt(ncol(dat) - 1))
mtry <- NULL

start_time <- Sys.time()
a <- SDTree(Y ~ ., dat, max_leaves = max_leaves, cp = cp, min_sample = min_sample, Q = Q, mtry = mtry, multicore = F)
end_time <- Sys.time()
end_time - start_time


ind <- lapply(1:nTree, function(x)sample(1:n, n, replace = T))
  #suppressWarnings({
  # estimating all the trees

X <- dat[, - ncol(dat)]
Y <- dat[, ncol(dat)]
data_list <- lapply(ind, function(i){
    return(list(X = X[i, ], Y = Y[i]))
})
res <- lapply(data_list, function(i)SDTree(x = i$X, y = i$Y, max_leaves = max_leaves, cp = cp, 
                                           min_sample = min_sample, Q = Q, mtry = mtry, multicore = F))


# mlbnenchmarks
library(microbenchmark)

data <- simulate_data_nonlinear(1, 200, 500, 5)
dat <- data.frame(X = data$X, Y = data$Y)
dat <- scale(dat, scale = T)
dat <- data.frame(dat)

microbenchmark(
  seq = SDTree(Y ~ ., data = dat, Q_type = 'DDL_trim', multicore = F),
  multicore = SDTree(Y ~ ., data = dat, Q_type = 'DDL_trim', multicore = T),
  times = 2
)

for(i in 1:100) a <- lapply(1:100, function(x) lapply(1:100, function(i) mean(sort(rnorm(10000)))))
a

res <- SDForest(Ozone ~ ., data, Q_type = 'DDL_trim', multicore = T)
res$predictions


plot(res$predictions, data$Ozone)


f <- function(x, a) {
  return(a * x)
}


ind <- 1:10
a <- 2
ind <- lapply(1:2, function(x)sample(1:n, n, replace = T))
f <- SDTree
X <- data[, -1]
Y <- data[, 1]
max_leaves <- 10
cp <- 0.01
min_sample <- 5
Q <- get_Q(X, 'pca', confounding_dim = 1)
mtry <- 2


    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, c("SDTree", "get_Q", "n_cores", "X", "Y", "max_leaves", 
                                  "cp", "min_sample", "Q", "mtry", "data.handler", "get_all_splitt", 
                                  "evaluate_splitt", "find_s", "predict_outsample", "traverse_tree", 
                                  "leave_names", "loss", "splitt_names"))
    res <- parLapply(cl = cl, X = ind, fun = function(i)SDTree(x = X[i, ], y = Y[i], max_leaves = max_leaves, cp = cp, min_sample = min_sample, 
                Q = Q, mtry = mtry, multicore = FALSE))
    parallel::stopCluster(cl = cl)
    print(res)


cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
parallel::clusterExport(cl, c("SDTree", "get_Q", "n_cores", "X", "Y", "max_leaves", 
                              "cp", "min_sample", "Q", "mtry", "data.handler", "get_all_splitt", 
                              "evaluate_splitt", "find_s", "predict_outsample", "traverse_tree", 
                              "leave_names"))
res <- parLapply(cl = cl, X = ind, fun = function(i)SDTree(x = X[i, ], y = Y[i], max_leaves = max_leaves, cp = cp, min_sample = min_sample, 
             Q = Q, mtry = mtry, multicore = FALSE))
parallel::stopCluster(cl = cl)
res

res =2

cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
res <- foreach(i = ind) %dopar% {
      fit <- SDTree(x = X[i, ], y = Y[i], max_leaves = max_leaves, cp = cp, min_sample = min_sample, 
             Q = Q, mtry = mtry, multicore = FALSE)
      return(fit)
    
    }
parallel::stopCluster(cl = cl)
res

res1 <- SDForest(Ozone ~ ., data, Q_type = 'DDL_trim', multicore = F, nTree = 2)
res2 <- SDForest(Ozone ~ ., data, Q_type = 'DDL_trim', multicore = T, nTree = 2)


source("utils.r")
data <- simulate_data_nonlinear(1, 10, 500, 5)


dat <- data.frame(X = data$X, Y = data$Y)

dat
Q <- get_Q(data$X, type = 'trim')
plot(rowSums(Q))

Q <- get_Q(scale(data$X), type = 'trim')
plot(rowSums(Q))


dat <- scale(dat, scale = T)

dat <- data.frame(dat)
dat


a <- SDTree(Y ~ ., dat, Q_type = 'DDL_trim', multicore = F, cp = 0.1)


res1 <- ranger(Y ~ ., data = dat, num.trees = 100)
res2 <- SDForest(Y ~ ., data = dat, nTree = 100, Q_type = 'DDL_trim', multicore = T, cp = 0.1)





dat$ranger <- res1$predictions
dat$SDForest <- res2$predictions

res2
res1





res1 <- ranger(Ozone ~ ., data = data, num.trees = 500)
res2 <- SDForest(Ozone ~ ., data = data, nTree = 500, Q_type = 'DDL_trim', multicore = T)

res2
res1

plot(dat)

dev.off()
plot(dat)

X <- as.matrix(dat[, -11])
plot(svd(X)$d)
Q <- get_Q(X, 'DDL_trim')

plot(rowSums(Q))

hist(Q)
heatmap(Q)
Q


plot(svd(Q)$d)


plot(svd(X)$d)
points(svd(Q %*% X)$d, col = 'red')

svd(Q %*% X)$d


E_tilde <- matrix(rowSums(Q))
Y_tilde <- Q %*% dat$Y

plot(E_tilde, Y_tilde)

c_hat <- fastLmPure(E_tilde, Y_tilde)$coefficients
c_hat

SDTree(Y ~ . , dat, Q_type = 'DDL_trim')

rnorm(100, 0, 100)


library(ggplot2)

n <- 100
p <- 200
q <- 5

sigma_nu <- 1
# strength of confounding
sigma_gamma <- sqrt(0.022)

# simulate data
set.seed(2023)
H <- matrix(rnorm(n * q), nrow = n)
Gamma <- matrix(rnorm(q * p, 0, sigma_gamma), nrow = q)
X <- matrix(rnorm(n * p, 0, sigma_nu), nrow = n)

X <- X + H %*% Gamma

X <- data$X

X <- scale(X)


decomp_conf <- svd(X)

Q <- get_Q(X, 'pca', confounding_dim = 1)
X_tilde <- Q %*% X


plot(decomp_conf$d)
points(c(0, svd(X_tilde)$d), col = 'red')






plot(svd(X)$d)
plot(svd(Q %*% X)$d)

