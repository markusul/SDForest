library(graphics)
library(ranger)
library(xgboost)

str(airquality)
data <- na.omit(airquality)
data$Ozone <- log(data$Ozone)
data <- scale(data)
data <- data.frame(data)
data

data$Day <- as.factor(data$Day)
str(data)
source("R/SDForest.r")

?get_Q

Q <- get_Q(data[1:80, -1], type = 'trim')

svd(Q)$d

Q_2 <- t(Q) %*% Q

sd_objective <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")

    grad <- - Q_2 %*% (labels - preds)
    hess <- Q_2
}

tree <- rpart(Ozone ~ ., data = data[1:80, ], control = rpart.control(xval = 20, cp = 0, minsplit = 5, minbucket = 2))
cptable <- tree$cptable
cptable


a = SDTree(Ozone ~ ., data = data[1:80, ], cp = 0.016)
a

a$predictions - predict(a, data[1:80,c(4, 2)])

res = cv.SDTree(Ozone ~ ., data = data[1:80, ], n_cv = 2)
res

res <- ranger(Ozone ~ ., data = data[1:80, ], num.trees = 100)
res

print(a)
plot(a)

b <- SDForest(Ozone ~ ., data = data[1:80, ])
b
b$f_X_hat
dev.off()
plot(data[1:80, 1])
points(a$predictions, col = 'red')
points(b$predictions, col = 'blue')

b$predictions - predict(b, data[1:80,])

predict(b, data[1:80,4])

!all(b$var.names %in% names(data))

library(rpart)
b = rpart(Ozone ~ ., data = data[1:80, ], control = rpart.control(xval = 20, cp = 0, minsplit = 2, minbucket = 2))
b$cptable

plot(diff(log(seq(1, 2, 0.1))))


print(a)
plot(a)

plot(data[1:80, 1])
points(a$f_X_hat, col = 'red')
points(b$f_X_hat, col = 'blue')

all(a$f_X_hat == b$f_X_hat)
print(a$tree, 'value', 's', 'j')
print(b$tree, 'value', 's', 'j')

print(a$vis_tree, 'value', 's', 'j', 'decision', 'label')
print(b$vis_tree, 'value', 's', 'j')

print(a)
plot(a)

SetEdgeStyle(a$tree, label = function(x) {x$decision})
SetNodeStyle(a$tree, label = function(x) {x$s})
plot(a$tree)


tree$Do(leave_names, filterFun = isLeaf)

tree[[3]]
print(tree, 'value', 's', 'j')
print(fit2$tree, 'value', 's', 'j')


fit1 <- SDTree(Y = data[1:80, 1], x = data[1:80, -1],min_sample = 2, cp = 0, multicore = F)
fit2 <- SDTree(y = data[1:80, 1], x = data[1:80, -1],cp = 0)
fit1$f_X_hat
fit2$f_X_hat

print(fit2)
plot(fit2)
predict(fit2, data[3,])
plot(data$Ozone)

fit2$tree[[1]]




tree_labels <- get_labels(fit2$vis_tree)
fit2$vis_tree
x = fit2
plot(x$vis_tree, layout = layout.reingold.tilford(x$vis_tree, root = 1), vertex.shape = tree_labels$node_shape, 
    vertex.label.color = 'black', vertex.size = 25, vertex.color = 'lightblue', edge.arrow.size = 0.5, 
    edge.label = tree_labels$arrow_lable, edge.label.color = 'black', edge.label.cex = 0.8)

plot(x$vis_tree)

fit1$f_X_hat - fit2$f_X_hat


print(fit$tree, 'value', 's', 'j')
data <- scale(data)


fit <- SDForest(Y = data[1:80, 1], X = data[1:80, -1], m = 100, nTree = 500, min_sample = 2, cp = 0, mtry = 4, deconfounding = T)
fit$f_X_hat
fit2 <- ranger(Ozone ~ ., data = data[1:80, ], num.trees = 100)

dtrain <- xgb.DMatrix(data = as.matrix(data[1:80, -1]), 
                      label = as.matrix(data[1:80, 1]))
fit3 <- xgb.train(data = dtrain, nrounds = 100)

pred <- predict(fit, data[81:111, -1])
pred2 <- predict(fit2, data[81:111, ])
pred3 <- predict(fit3, xgb.DMatrix(data = as.matrix(data[81:111, -1])))

mean((data[81:111, 1] - pred)**2)
mean((data[81:111, 1] - pred2$predictions)**2)
mean((data[81:111, 1] - pred3)**2)

Q_test <- get_Q(data[81:111, -1], type = 1)
mean((Q_test %*% data[81:111, 1] - Q_test %*% pred)**2)
mean((Q_test %*% data[81:111, 1] - Q_test %*% pred2$predictions)**2)
mean((Q_test %*% data[81:111, 1] - Q_test %*%pred3)**2)


?ranger

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

simulate_data_nonlinear_E <- function(q, p, n, m){
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

    b <- rbinom(n, 1, 0.5)
    e <- rnorm(n, -2, 1)
    e[b == 1] <- rnorm(sum(b), 2, 1)
    E[, 1] <- e

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

data <- simulate_data_nonlinear_E(1, 20, 100, 5)
data

plot(svd(scale(data$X))$d)


plot(density(rnorm(100, 0, 1)))
b <- rbinom(100, 1, 0.5)
e <- rnorm(100, -2, 1)
e[b == 1] <- rnorm(sum(b), 2, 1)
e

plot(density(rnorm(100, -0.5, 1) + rnorm(100, 0.5, 1)))
plot(density(e))



simulate_data_nonlinear_D <- function(q, p, n, m, d = 1){
    #simulate data with confounding and non-linear f_X
    # q: number of confounding covariates in H
    # p: number of covariates in X
    # n: number of observations
    # m: number of covariates with a causal effect on Y

    # random confounding covariates H
    H <- matrix(rnorm(n * q, 0, 1), nrow = n)

    # random correlation matrix cov(X, H)
    Gamma <- matrix(rnorm(q * p, 0, 1), nrow = q)
    non_confounded_covs <- sample(1:p, floor((1-d) * p))
    print(non_confounded_covs)
    Gamma[, non_confounded_covs] <- 0


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



data <- simulate_data_nonlinear_D(1, 20, 100, 5, d = 0.1)
data

quantile(sv$d, 0.5)

get_Q <- function(X, type, q_hat = 0){
    # X: covariates
    # type: type of deconfounding
    modes <- c('trim' = 1, 'DDL_trim' = 2, 'pca' = 3, 'no_deconfounding' = 4)
    if(!(type %in% names(modes))) stop(paste("type must be one of:", paste(names(modes), collapse = ', ')))

    # number of observations
    n <- dim(X)[1]

    # calculate deconfounding matrix
    sv <- svd(X)
    tau <- median(sv$d)
    D_tilde <- unlist(lapply(sv$d, FUN = function(x)min(x, tau))) / sv$d

    Q <- switch(modes[type], sv$u %*% diag(D_tilde) %*% t(sv$u), # trim
                            diag(n) - sv$u %*% diag(1 - D_tilde) %*% t(sv$u), # DDL_trim
                            { # pca
                                d_pca <- sv$d
                                if(q_hat <= 0) stop("the assumed confounding dimension q_hat must be larger than zero")
                                d_pca[1:q_hat] <- 0
                                print('aa')
                                sv$u %*% diag(d_pca) %*% t(sv$u)
                            },
                            diag(n)) # no_deconfounding
    return(Q)
}

print(paste(c(1, 2, 3, 4)))


X <- matrix(rnorm(100 * 20), nrow = 100)
get_Q(X, 'pa')

paste(names(modes), collapse = ', ')

sv <- svd(X)

type <- 'pc'
q_hat <- 0
!(type %in% names(modes))


modes <- c('trim' = 1, 'DDL_trim' = 2, 'pca' = 3, 'no_deconfounding' = 4)

Q <- switch(modes[type], sv$u %*% diag(D_tilde) %*% t(sv$u), 
                         diag(n) - sv$u %*% diag(1 - D_tilde) %*% t(sv$u), 
                         {
                            d_pca <- sv$d
                            if(q_hat <= 0) stop("the assumed confounding dimension q_hat must be larger than zero")
                            d_pca[1:q_hat] <- 0
                            print('aa')
                            sv$u %*% diag(d_pca) %*% t(sv$u)
                         },
                         diag(n))
