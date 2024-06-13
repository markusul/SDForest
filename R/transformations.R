#' Estimation of spectral transformation
#' 
#' Estimates the spectral transformation Q for spectral deconfounding by 
#' shrinking the leading singular values of the covariates.
#' @importFrom Rdpack reprompt
#' @references
#'   \insertAllCited{}
#' @author Markus Ulmer
#' @param X Numerical covariates of class \code{matrix}.
#' @param type Type of deconfounding, one of 'trim', 'pca', 'no_deconfounding'. 
#' 'trim' corresponds to the Trim transform \insertCite{Cevid2020SpectralModels}{SDForest} 
#' as implemented in the Doubly debiased lasso \insertCite{Guo2022DoublyConfounding}{SDForest}, 
#' 'pca' to the PCA transformation\insertCite{Paul2008PreconditioningProblems}{SDForest} 
#' and 'no_deconfounding' to the Identity.
#' @param trim_quantile Quantile for Trim transform, only needed for trim.
#' @param q_hat Assumed confounding dimension, only needed for pca.
#' @param gpu If \code{TRUE}, the calculations are performed on the GPU. 
#' If it is properly set up.
#' @return Q of class \code{matrix}, the spectral transformation matrix.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(50 * 20), nrow = 50)
#' Q_trim <- get_Q(X, 'trim')
#' Q_pca <- get_Q(X, 'pca', q_hat = 5)
#' Q_plain <- get_Q(X, 'no_deconfounding')
#' @export
get_Q <- function(X, type, trim_quantile = 0.5, q_hat = 0, gpu = FALSE){
  if(gpu) ifelse(GPUmatrix::installTorch(), 
                 gpu_type <- 'torch', 
                 gpu_type <- 'tensorflow')
  
  if(type == 'no_deconfounding') {
    Q <- diag(nrow(X))
    if(gpu) Q <- gpu.matrix(Q, type = gpu_type)
    return(Q)
  }
  
  
  svd_error <- function(X, f = 1, count = 1){
    tryCatch({
      svd(X * f)
    }, error = function(e) {
      warning(paste(e, ':X multipied by number close to 1'))
      if(count > 5) stop('svd did not converge')
      return(svd_error(X, 1 + 0.0000000000000001 * 10 ^ count, count + 1))})
  }
  
  if(ncol(X) == 1){
    warning('only one covariate, no deconfounding possible')
    Q <- diag(nrow(X))
    if(gpu) Q <- gpu.matrix(Q, type = gpu_type)
    return(Q)
  }
  
  modes <- c('trim' = 1, 'pca' = 2, 'no_deconfounding' = 3)
  if(!(type %in% names(modes))) stop(paste("type must be one of:", 
                                           paste(names(modes), collapse = ', ')))
  
  # number of observations
  n <- dim(X)[1]
  
  # calculate deconfounding matrix
  sv <- svd_error(X)
  
  tau <- quantile(sv$d, trim_quantile)
  D_tilde <- unlist(lapply(sv$d, FUN = function(x)min(x, tau))) / sv$d
  D_tilde[is.na(D_tilde)] <- 1
  
  U <- sv$u
  if(gpu) U <- gpu.matrix(U, type = gpu_type)
  
  switch(modes[type], 
         diag(n) - U %*% diag(1 - D_tilde) %*% t(U), # DDL_trim
         { # pca
           d_pca <- rep(1, length(sv$d))
           if(q_hat <= 0) 
             stop("the assumed confounding dimension must be larger than zero")
           d_pca[1:q_hat] <- 0
           diag(n) - U %*% diag(1 - d_pca) %*% t(U)
            },
         diag(n) # no_deconfounding
         )
}

#' Estimation of anchor transformation
#' 
#' Estimates the anchor transformation for the Anchor-Objective. 
#' The anchor transformation is \eqn{W = I-(1-\sqrt{\gamma}))\Pi_A}, 
#' where \eqn{\Pi_A = A(A^TA)^{-1}A^T}. For \eqn{\gamma = 1} this is just the identity. 
#' For \eqn{\gamma = 0} this corresponds to residuals after orthogonal projecting onto A.
#' For large \eqn{\gamma} this is close to the orthogonal projection onto A, scaled by \eqn{\gamma}.
#' The estimator \eqn{\text{argmin}_f ||W(Y - f(X))||^2} corresponds to the Anchor-Regression Estimator 
#' \insertCite{Rothenhausler2021AnchorCausality}{SDForest}, \insertCite{Buhlmann2020InvarianceRobustness}{SDForest}.
#' @importFrom Rdpack reprompt
#' @references
#'   \insertAllCited{}
#' @author Markus Ulmer
#' @param A Numerical Anchor of class \code{matrix}.
#' @param gamma Strength of distributional robustness, \eqn{\gamma \in [0, \infty]}.
#' @param intercept Logical, whether to include an intercept in the anchor.
#' @param gpu If \code{TRUE}, the calculations are performed on the GPU. 
#' If it is properly set up.
#' @return W of class \code{matrix}, the anchor transformation matrix.
#' @examples
#' set.seed(1)
#' n <- 50
#' X <- matrix(rnorm(n * 1), nrow = n)
#' Y <- 3 * X + rnorm(n)
#' W <- get_W(X, gamma = 0)
#' resid <- W %*% Y
#' @export
get_W <- function(A, gamma, intercept = FALSE, gpu = FALSE){
  if(intercept) A <- cbind(1, A)
  if(ncol(A) > nrow(A)) stop('A must have full rank!')
  if(gamma < 0) stop('gamma must be non-negative')
  
  if(gpu) ifelse(GPUmatrix::installTorch(), 
                 gpu_type <- 'torch', 
                 gpu_type <- 'tensorflow')  
  if(gpu) A <- gpu.matrix(A, type = gpu_type)
  
  Q_prime <- qr.Q(qr(A))
  Pi_A <- tcrossprod(Q_prime)
  diag(nrow(A)) - (1-sqrt(gamma)) * Pi_A
}