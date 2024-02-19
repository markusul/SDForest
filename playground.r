source("R/SDForest.r")
library(ggplot2)

p <- 20
n <- 20

set.seed(2024)
data <- simulate_data_nonlinear(20, p, n, 4)





# complexity of f_X
complexity <- 5
# random parameter for fourier basis
beta <- runif(m * complexity * 2, -1, 1)
# random confounding covariates H
H <- matrix(rnorm(n * q, 0, 1), nrow = n)
# random correlation matrix cov(X, H)
Gamma <- matrix(rnorm(q * p, 0, 1), nrow = q)
if(!is.null(eff)){
    non_effected <- p - eff
    if(non_effected <= 0) stop('eff must be smaller than p or NULL'){
      Gamma[, sample(1:p, non_effected)] <- 0
    }
}
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

