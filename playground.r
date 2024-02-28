source("R/SDForest.r")
library(ggplot2)

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
