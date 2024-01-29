
simulate_data_step <- function(q, p, n, m, make_tree = F){
    # q: number of confounding covariates in H
    # p: number of covariates in X
    # n: number of observations
    # m: number of partitions

    # minimum number of observations for split
    min_sample <- 2

    # random confounding covariates H
    H <- matrix(rnorm(n * q), nrow = n)

    # random correlation matrix cov(X, H)
    Gamma <- matrix(rnorm(q * p), nrow = q)

    # random coefficient vector delta
    delta <- rnorm(q, 0, 10)

    # relevant covariates
    js <- c()

    # generate X
    if(q == 0){
        # no confounding covariates
        X <- matrix(rnorm(n * p), nrow = n)
    }else{
        # confounding covariates
        X <- H %*% Gamma + matrix(rnorm(n * p), nrow = n)
    }

    # generate tree
    if(make_tree){
        tree <- Node$new(name = '1', value = 0)
    }

    # partitions of observations
    index <- list(1:n)
    for (i in 1:m){
        # get number of observations in each partition
        samples_per_part <- unlist(lapply(index, function(x)length(x)))

        # get potential splits (partitions with enough observations)
        potential_splitts <- which(samples_per_part >= min_sample)

        # sample potential partition to split
        # probability of each partition is proportional to number of observations
        if(length(potential_splitts) == 1){
            branch <- potential_splitts
        }else {
            branch <- sample(potential_splitts, 1,
                            prob = samples_per_part[potential_splitts]/sum(samples_per_part[potential_splitts]))
        }

        # sample covariate to split on
        j <- sample(1:p, 1)
        js <- c(j, js)

        # sample split point
        potential_s <- X[index[[branch]], j]
        s <- rbeta(1, 2, 2) * (max(potential_s) - min(potential_s)) + min(potential_s)

        # split partition
        index <- append(index, list(index[[branch]][X[index[[branch]], j] > s]))
        index[[branch]] <- index[[branch]][X[index[[branch]], j] <= s]

        # add split to tree
        if(make_tree){
            if(tree$height == 1){
                leave <- tree
            }else{
                leaves <- tree$leaves
                leave <- leaves[[which(tree$Get('name', filterFun = isLeaf) == branch)]]
            }
            leave$j <- j
            leave$s <- s
            leave$AddChild(branch, value = 0)
            leave$AddChild(i + 1, value = 0)
        }
     }

    # sample means per partition
    f_X_means <- runif(length(index), -50, 50)

    # generate f_X
    f_X <- rep(0, n)
    for (i in 1:length(index)){
        f_X[index[[i]]] <- f_X_means[i]
    }

    # generate Y
    if(q == 0){
        # no confounding covariates
        Y <- f_X + rnorm(n)
    }else{
        # confounding covariates
        Y <- f_X + H %*% delta + rnorm(n)
    }

    # return data
    if(make_tree){
        # add leave values to tree
        for(l in tree$leaves){
            l$value <- f_X_means[as.numeric(l$name)]
        }
        return(list(X = X, Y = Y, f_X = f_X, j = js, tree = tree, index = index))
    }

    # return data
    return(list(X = X, Y = Y, f_X = f_X, j = js, index = index))
}

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

source("R/SDForest.r")

#source('utils.r')
p <- 50
n <- 200
data <- simulate_data_nonlinear(10, p, n, 1)


dat <- data.frame(X = data$X, Y = data$Y)
#dat <- scale(dat, scale = T)
dat <- data.frame(dat)

a <- SDTree(Y ~ ., dat, Q_type = 'DDL_trim', cp = 0.01, max_leaves = 400)
plot(a$var_imp)


b <- SDForest(Y ~ ., dat, Q_type = 'DDL_trim', cp = 0.01, max_leaves = 400, nTree = 100)
f <- condDependence(b, dat, data$j[1])
plot(b$var_importance)


plot(data$X[, data$j[1]], data$Y)
points(data$X[, data$j[1]], a$predictions, col = 'red', pch = 20)
points(data$X[, data$j[1]], b$predictions, col = 'blue', pch = 20)
points(data$X[, data$j[1]], data$f_X, col = 'green', pch = 20)

library(ranger)
c <- ranger(Y ~ ., data = dat, num.trees = 100, importance = 'impurity')
points(data$X[, data$j[1]], c$predictions, col = 'purple', pch = 2, cex = 0.5)

plot(c$variable.importance/max(c$variable.importance), ylim = c(0, 1))
points(b$var_imp/max(b$var_imp), col = 'red', pch = 20)
points(a$var_imp/max(a$var_imp), col = 'blue', pch = 20)
points(1, data$j[1], col = 'green', pch = 20)


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

