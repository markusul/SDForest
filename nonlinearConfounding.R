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