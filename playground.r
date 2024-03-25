source("R/SDForest_gpu.r")

library('ranger')

m <- 5
p <- 1000
n <- 1000
q <- 4

data <- simulate_data_nonlinear(q, p, n, m)

X <- data$X
Y <- data$Y


a <- Sys.time()
fit1 <- SDTree(x = X, y = Y)
b <- Sys.time()

source("R/SDForest.r")
c <- Sys.time()
fit2 <- SDTree(x = X, y = Y)
d <- Sys.time()

b - a
d - c

max(fit1$predictions - fit2$predictions)





lm.fit(X_gpu, Y_gpu)$coefficients
dat <- data.frame(X, Y)
GPUglm(dat, formula = Y ~ ., family = gaussian())$coefficients

#linear model using glm.fit.gpu
utils::data(anorexia, package = "MASS")
anorex_glm <- glm(Postwt ~ Prewt + Treat + offset(Prewt),
                  family = gaussian(), data = anorexia)
summary(anorex_glm)
x <- model.matrix(~Treat+Prewt,data=anorexia)
y <- as.matrix(anorexia$Postwt)
s1_glm <- glm.fit(x=x,y=y)
s1_gpu <- glm.fit.GPU(x=x,y=y)

a <- Sys.time()
fit2 <- SDTree(x = data$X, y = data$Y)
b <- Sys.time()

b -a

fit2$predictions == fit$predictions


install.packages("tensorflow")
library(tensorflow)
install_tensorflow(version = "nightly-gpu")
library(GPUmatrix)

X_gpu <- gpu.matrix(X, device = 'cuda')

svd(X_gpu)

?gpu.matrix

X_gpu@gm$device

svd(X_gpu)


X_gpu %*% t(X_gpu)


Gm <- gpu.matrix(c(1:20)+40,10,2, device = 'cuda')

m <- matrix(rnorm(5000), 1000)
mg <- gpu.matrix(m)
Gm
mg

svd(mg)


m <- matrix(c(1:20)+40,10,2)
Gm <- gpu.matrix(c(1:20)+40,10,2)

head(tcrossprod(m),1)
head(tcrossprod(Gm),1)



install.packages("torch")
library(torch)
install_torch()

library(GPUmatrix)

n <- 1000
p <- 1000
X <- matrix(rnorm(n * p), nrow = n)
b <- matrix(rnorm(p), nrow = p)
G <- gpu.matrix(X, device= 'cuda')
g <- gpu.matrix(b, device= 'cuda')

a <- Sys.time()
x <- lapply(1:100, function(rep) X - b %*% (t(b) %*% X))
b <- Sys.time()
x <- lapply(1:100, function(rep) G - g %*% (t(g) %*% G))
c <- Sys.time()

b - a
c - b


get_all_splitt <- function(branch, X, n, min_sample, p, E){
  # finds the best splitts for every covariate in branch
  # returns the best splitpoint for every covariate and the resulting loss decrease

  # observations belonging to branch
  index <- which(E[, branch] == 1)

  X_branch <- X[index, ]

  # all possible split points
  s <- find_s(X_branch, min_sample, p)

  

  res <- lapply(1:p, function(j) {lapply(s[, j], function(x) {
            X_branch_j <- if(p == 1) X_branch else X_branch[, j]
            eval <- evaluate_splitt(branch = branch, j = j, 
              s = x, index = index, X_branch_j = X_branch_j, n = n, 
              min_sample = min_sample)
            return(eval)})})

  res <- lapply(res, function(x)do.call(rbind, x))
  res_opt <- lapply(res, function(x)x[which.max(x[, 1]), ])
  res_opt <- do.call(rbind, res_opt)

  return(matrix(unlist(res_opt), ncol = 4, byrow = F))
}

evaluate_splitt <- function(branch, j, s, index, X_branch_j, n, min_sample){
  # evaluate a split at partition branch on covariate j at the splitpoint s
  # index: index of observations in branch
  # dividing observation in branch
  if (sum(X_branch_j > s) < min_sample | sum(X_branch_j <= s) < min_sample){
    return(list('dloss' = 0, j = j, s = s, branch = branch))
  }

  e_next <- e_next * 0
  e_next[index[X_branch_j > s]] <- 1

  u_next_prime <- Q_temp %*% e_next
  u_next_size <- sum(u_next_prime ** 2)

  dloss <- crossprod(u_next_prime, Y_tilde)**2 / u_next_size

  return(list('dloss' = as.numeric(dloss), j = j, s = s, branch = branch))
}