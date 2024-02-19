# nonlinear confounding
library(splines)

n <- 200
p <- 600

E <- matrix(rnorm(n*p), ncol=p)

H <- rnorm(n)

# size of the b-spline basis.
df <- 15

#use b-spline basis to generate "random function" from R to R^p
bb <- bs(H, df = df, intercept=F)

# visualize a random function generated in this way.
w <- runif(df)
b <- bb%*%w
plot(H, b)

Gam <- matrix(runif(p*df, min=-1, max=1), ncol=df)

X <- bb%*%t(Gam)+E

svd(X)$d
plot(svd(X)$d)

#we see many spikes even though q=1


library(EQL)


f_hermite <- function(x, beta, js){
    # function to generate f_X
    # x: covariates
    # beta: parameter vector
    # js: relevant covariates

    # number of relevant covariates
    m <- length(js)

    # complexity of f_X
    complexity <- length(beta) / m

    if(is.null(dim(beta))) beta <- matrix(beta)

    # calculate f_X
    do.call(sum, lapply(1:m, function(i) {
        j <- js[i]
        # select beta for covariate j
        res <- lapply(1:complexity, function(k) {
            beta[k, i] * hermite(x[j], k-1)
        })
        Reduce('+', res)
    }))
}

p <- 500
q <- 1
n <- 500
H <- matrix(rnorm(n*q), ncol=q)
df <- 5

Betas <- replicate(p, matrix(runif(q*df, -1, 1), nrow = df), simplify = F)

X <- sapply(Betas, function(beta) {apply(H, 1, function(h) f_hermite(h, beta, 1:q))})
X <- X + matrix(rnorm(n*p), ncol = p)

plot(svd(X)$d)

