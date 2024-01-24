library(graphics)
library(ranger)
library(xgboost)
library(rpart)

str(airquality)
data <- na.omit(airquality)
data$Ozone <- log(data$Ozone)
data <- scale(data)
data <- data.frame(data)
data



source("R/SDForest.r")

source('utils.r')
data <- simulate_data_nonlinear(1, 200, 300, 5)

dat <- data.frame(X = data$X, Y = data$Y)
dat <- scale(dat, scale = T)
dat <- data.frame(dat)

start_time <- Sys.time()
a <- SDTree(Y ~ ., dat, Q_type = 'DDL_trim', multicore = F, cp = 0.01, max_leaves = 300, mtry = 180)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
a <- SDForest(Y ~ ., dat, Q_type = 'DDL_trim', multicore = T, mtry = 180, cp = 0.01, max_leaves = 300, nTree = 100)
end_time <- Sys.time()
end_time - start_time



b <- a$predictions

plot(b, a$predictions)


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
