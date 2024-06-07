#' Simulate data with linear confounding and non-linear causal effect
#' 
#' Simulation of data from a confounded non-linear model. 
#' The data generating process is given by:
#' \deqn{Y = f(X) + H \delta + \nu}
#' \deqn{X = H \Gamma + E}
#' where \eqn{f(X)} is a random function on the fourier basis
#' with a subset of size m covariates \eqn{X_j} having a causal effect on \eqn{Y}.
#' \deqn{f(x_i) = \sum_{j = 1}^p 1_{j \in js} \sum_{k = 1}^K \beta_{j, k, 1}^{(1)} \cos(0.1 k x_j) + 
#' \beta_{j, k, 2}^{(2)} \sin(0.1 k x_j)}
#' \eqn{E}, \eqn{\nu} are random error terms and 
#' \eqn{H \in \mathbb{R}^{n \times q}} is a matrix of random confounding covariates.
#' \eqn{\Gamma \in \mathbb{R}^{q \times p}} and \eqn{\delta \in \mathbb{R}^{q}} are random coefficient vectors.
#' For the simulation, all the above parameters are drawn from a standard normal distribution, except for 
#' \eqn{\nu} which is drawn from a normal distribution with standard deviation 0.1.
#' The parameters \eqn{\beta} are drawn from a uniform distribution between -1 and 1.
#' @author Markus Ulmer
#' @param q number of confounding covariates in H
#' @param p number of covariates in X
#' @param n number of observations
#' @param m number of covariates with a causal effect on Y
#' @param eff the number of affected covariates in X by the confounding, if NULL all covariates are affected
#' @return a list containing the simulated data:
#' \item{X}{a matrix of covariates}
#' \item{Y}{a vector of responses}
#' \item{f_X}{a vector of the true function f(X)}
#' \item{j}{the indices of the causal covariates in X}
#' \item{beta}{the parameter vector for the function f(X)}, see \code{\link{f_four}}
#' \item{H}{the matrix of confounding covariates}
#' @seealso \code{\link{f_four}}
#' @export 
simulate_data_nonlinear <- function(q, p, n, m, eff = NULL){
  #simulate data with confounding and non-linear f_X
  # q: number of confounding covariates in H
  # p: number of covariates in X
  # n: number of observations
  # m: number of covariates with a causal effect on Y

  # complexity of f_X (number of fourier basis functions) K
  complexity <- 5
  # random parameter for fourier basis
  beta <- runif(m * complexity * 2, -1, 1)

  # random confounding covariates H
  H <- matrix(rnorm(n * q, 0, 1), nrow = n)

  # random correlation matrix cov(X, H)
  Gamma <- matrix(rnorm(q * p, 0, 1), nrow = q)

  if(!is.null(eff)){
    non_effected <- p - eff
    if(non_effected <= 0) stop('eff must be smaller than p or NULL')

    Gamma[, sample(1:p, non_effected)] <- 0
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
  
  #return data
  return(list(X = X, Y = Y, f_X = f_X, j = js, beta = beta, H = H))
}

#' Function of x on a fourier basis
#' 
#' Function of x on a fourier basis with a subset of covariates 
#' having a causal effect on Y using the parameters beta.
#' The function is given by:
#' \deqn{f(x_i) = \sum_{j = 1}^p 1_{j \in js} \sum_{k = 1}^K \beta_{j, k, 1}^{(1)} \cos(0.1 k x_j) +
#' \beta_{j, k, 2}^{(2)} \sin(0.1 k x_j)}
#' @author Markus Ulmer
#' @param x a vector of covariates
#' @param beta the parameter vector for the function f(X)
#' @param js the indices of the causal covariates in X
#' @return the value of the function f(x)
#' @seealso \code{\link{simulate_data_nonlinear}}
#' @export
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
    do.call(sum, lapply(1:complexity, function(k) 
      beta[beta_ind[1 + (k-1) *2]] * sin(k * 0.1 * x[j]) + 
        beta[beta_ind[2 + (k-1) *2]] * cos(k * 0.1 * x[j])))
  }))
}