source("R/SDForest.r")
library(ggplot2)

m <- 5
p <- 100
n <- 100
q <- 4


complexity <- 5
# random parameter for fourier basis
beta <- runif(m * complexity * 2, -1, 1)
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
# generate f_X
f_X <- apply(X, 1, function(x) f_four(x, beta, js))

# generate Y
Y <- f_X + H %*% delta + rnorm(n, 0, 0.1)


H_test <- matrix(rnorm(n * q, 0, 1), nrow = n)
H_test[, 1:2] <- H_test[, 1:2] + 10
E_test <- matrix(rnorm(n * p, 0, 1), nrow = n)
X_test <- H_test %*% Gamma + E_test
f_X_test <- apply(X_test, 1, function(x) f_four(x, beta, js))
Y_test <- f_X_test + H_test %*% delta + rnorm(n, 0, 0.1)

library(ranger)
train_data <- data.frame(X, Y)

res0 <- ranger(Y ~ ., train_data, num.trees = 100, importance = 'impurity', mtry = floor(0.9 * ncol(X)))
res1 <- SDForest(Y ~ ., train_data)
res2 <- SDForest(Y ~ ., train_data, A = H[, 1:2], gamma = 2, Q_type = 'no_deconfounding')
res3 <- SDForest(Y ~ ., train_data, A = H[, 1:2], gamma = 2, Q_type = 'trim')
res4 <- SDForest(Y ~ ., train_data, A = H[, 1:2], gamma = 0.5, Q_type = 'no_deconfounding')
res5 <- SDForest(Y ~ ., train_data, A = H[, 1:2], gamma = 0.5, Q_type = 'trim')

res0_test <- predict(res0, data.frame(X_test))$predictions
res1_test <- predict(res1, data.frame(X_test))
res2_test <- predict(res2, data.frame(X_test))
res3_test <- predict(res3, data.frame(X_test))
res4_test <- predict(res4, data.frame(X_test))
res5_test <- predict(res5, data.frame(X_test))

mean((res0$predictions - f_X)**2)
mean((res1$predictions - f_X)**2)
mean((res2$predictions - f_X)**2)
mean((res3$predictions - f_X)**2)
mean((res4$predictions - f_X)**2)
mean((res5$predictions - f_X)**2)

print('a')
mean((res0$predictions - Y)**2)
mean((res1$predictions - Y)**2)
mean((res2$predictions - Y)**2)
mean((res3$predictions - Y)**2)
mean((res4$predictions - Y)**2)
mean((res5$predictions - Y)**2)

print('b')
mean((res0_test - f_X_test)**2)
mean((res1_test - f_X_test)**2)
mean((res2_test - f_X_test)**2)
mean((res3_test - f_X_test)**2)
mean((res4_test - f_X_test)**2)
mean((res5_test - f_X_test)**2)

print('c')
mean((res0_test - Y_test)**2)
mean((res1_test - Y_test)**2)
mean((res2_test - Y_test)**2)
mean((res3_test - Y_test)**2)
mean((res4_test - Y_test)**2)
mean((res5_test - Y_test)**2)





set.seed(22)
n <- 12
p <- 2
X <- matrix(rnorm(n * p), ncol = p)

y <- X %*% rnorm(p)

X[2:n, 1] <- 1
svd(X)$d

get_Q(X, 'trim')

res <- SDForest(x = X, y = y)
res



p <- 5000
n <- 5000

set.seed(2024)
data <- simulate_data_nonlinear(20, p, n, 4)


start <- Sys.time()
svd(data$X)$d
end <- Sys.time()
end - start

fit <- SDTree(x = data$X, y = data$Y)


cv.SDTree(x = data$X, y = data$Y)
prune(fit, 0.52)

start <- Sys.time()
for(i in 1:10) fit <- SDForest(x = data$X, y = data$Y)
end <- Sys.time()
(end - start) / 10
# Time difference of 0.8825729 mins


A  <- matrix(rnorm(4*4), 4, 4)
A <- rbind(A, A[1, ])
A
unique(A)
find_s(A, 2, 4)


a <- regPath(fit, oob = T)
b <- stabilitySelection(fit)

multicore <- T

p <- 500
n <- 500
q <- 20

n_test <- 500

set.seed(2024)
data <- simulate_data_nonlinear(q, p, n + n_test, 4)
data_test <- data
data_test$Y <- data_test$Y[(n+1):(n+n_test)]
data_test$X <- data_test$X[(n+1):(n+n_test),]
data_test$f_X <- data_test$f_X[(n+1):(n+n_test)]

data$X <- data$X[1:n,]
data$Y <- matrix(data$Y[1:n])
data$f_X <- data$f_X[1:n]


  input_data <- data.handler(x = data$X, y = data$Y)
  X <- input_data$X
  Y <- input_data$Y

  # number of observations
  n <- nrow(X)
  # number of covariates
  p <- ncol(X)


  ind <- lapply(1:100, function(x)sample(1:n, n, replace = T))

for(i in ind) svd(X[i, ])$d

m <- cor(X[i, ])

diag(m) <- 0

sort(m, decreasing = T)[1:60]
svd(X[i, ]*1.000000000000001)$d

a <- tryCatch({
  svd(X[i, ])$d
}, error = function(e) {
    warning(paste(e, ':X multipied by number close to 1'))
    return(svd(X[i, ] * 1.000000000000001)$d)})
a

p <- 10
Q <- matrix(rnorm(p^2), p, p)
Q

E <- matrix(sample(c(0, 1), 10, replace = T))
E <- cbind(E, as.numeric(E == 0))
E

rowSums(Q[, as.logical(E[, 2])])

Q%*%E


p <- 1000
n <- 1000

set.seed(2024)
data <- simulate_data_nonlinear(20, p, n, 4)
Q <- get_Q(data$X, 'trim')


start <- Sys.time()
E_tilde <- Q %*% matrix(1, ncol(Q), 1)
end <- Sys.time()

(end - start) * 100



Q <- matrix(rnorm(4), 2)
norm(Q %*% c(1, 0), type = '2')

sqrt(sum(rnorm()**2))

n <- 10
u <-  matrix(rnorm(n), n, 1)
Q <- matrix(rnorm(n**2), n, n)

a <- Sys.time()
x <- SMUT::eigenMapMatMult(u, SMUT::eigenMapMatMult(t(u), Q))
b <- Sys.time()
x <- u %*% (t(u) %*% Q)
d <- Sys.time()
#x <- u %*% t(u) %*% Q
e <- Sys.time()

(b - a)
(d - b)
(e - d)

Q <- matrix(1, 10, 10)
Q
E_tilde <- rowSums(Q)

sum((E_tilde / sqrt(sum(E_tilde**2)))**2)



SMUT::eigenMapMatMult(u, SMUT::eigenMapMatMult(t(u), Q))


source("R/SDForest.r")
A <- matrix(rnorm(10 * 5), ncol = 5)
A

A

d <- qr(A)
d$qr

Q_prime <- qr.Q(d)

Pi_A <- tcrossprod(Q_prime)


W <- diag(10) - Pi_A

t(W) == W


library(matlib)
QR(A)

source("R/SDForest.r")
data <- simulate_data_nonlinear(5, 100, 100, 4)

X <- data$X
A <- data$H[, 1:2]
Y <- data$Y

fit <- SDForest(x = X, y = Y)
fit2 <- SDForest(x = X, y = Y, A = A)

mean((fit$predictions - data$f_X)**2)
mean((fit2$predictions - data$f_X)**2)

X <- data$X

plot(svd(X)$d)

Q_prime <- qr.Q(qr(A))
Pi_A <- tcrossprod(Q_prime)
W <- diag(nrow(A)) - (1-sqrt(gamma)) * Pi_A


X_prime <- W %*% X

plot(svd(X_prime)$d)

Q <- get_Q(X_prime, 'trim')

plot(svd(Q %*% W %*% X)$d)

svd(get_Q(X, 'trim') %*% X)$d



