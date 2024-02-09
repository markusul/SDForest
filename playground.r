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
        do.call(sum, lapply(1:complexity, function(k) beta[beta_ind[1 + (k-1) *2]] * sin(k * 0.2 * x[j]) + beta[beta_ind[2 + (k-1) *2]] * cos(k * 0.2 * x[j])))
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

source("R/SDForest.r")

p <- 100
n <- 100
data <- simulate_data_nonlinear(20, p, n, 4)


dat <- data.frame(X = data$X, Y = data$Y)
#dat <- scale(dat, scale = T)
dat <- data.frame(dat)

start_time <- Sys.time()
b <- SDForest(Y ~ ., dat,min_sample = 1)
print(Sys.time() - start_time)


plot(b$predictions)
points(dat$Y, col = 'red', pch = 20)
points(data$f_X, col = 'green', pch = 20)
b

f <- condDependence(b, dat, data$j[1])
plot(f)

plot(a$var_importance)

true_f <- true_function(data$beta, data$j)
predict(true_f, data$X)

dep_f <- condDependence(true_f, data$X, data$j[1])
plot(dep_f)

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

