#' COmputes retval
#' 
#' @param x Vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#' @param n  Number of observations.
#' @param mean Mean vector, default is rep(0, length = ncol(x)).
#' @param sigma Covariance matrix, default is diag(ncol(x)).
#' @param log	Logical; if TRUE, densities d are given as log(d).
#' @param trustme	Logical; skip checks for full speed. You're on your own here!
#' @param method	Matrix decomposition used to determine the matrix root of sigma, possible methods are eigenvalue decomposition ("eigen", default), singular value decomposition ("svd"), and Cholesky decomposition ("chol").
#' @param pre0.9_9994	# Logical; if FALSE, the output produced in mvtnorm versions up to 0.9-9993 is reproduced. In 0.9-9994, the output is organized such that rmvnorm(10,...) has the same first ten rows as rmvnorm(100, ...) when called with the same seed.
#' @return retval
#' @keywords internal
rmvnorm_withGivenLatent <-
function (n, latent, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
          method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE) 
{
#   x Vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#   
#   n	Number of observations.
#   
#   mean Mean vector, default is rep(0, length = ncol(x)).
#   
#   sigma Covariance matrix, default is diag(ncol(x)).
#   
#   log	Logical; if TRUE, densities d are given as log(d).
#   
#   trustme	Logical; skip checks for full speed. You're on your own here!
# 
#   method	Matrix decomposition used to determine the matrix root of sigma, possible methods are eigenvalue decomposition ("eigen", default), singular value decomposition ("svd"), and Cholesky decomposition ("chol").
# 
#   pre0.9_9994	# Logical; if FALSE, the output produced in mvtnorm versions up to 0.9-9993 is reproduced. In 0.9-9994, the output is organized such that rmvnorm(10,...) has the same first ten rows as rmvnorm(100, ...) when called with the same seed.
  
  
  
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                   check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != nrow(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  sigma1 <- sigma
  dimnames(sigma1) <- NULL
  if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
    warning("sigma is numerically not symmetric")
  }
  method <- match.arg(method)
  if (method == "eigen") {
    ev <- eigen(sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigma is numerically not positive definite")
    }
    retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
  }
  else if (method == "svd") {
    sigsvd <- svd(sigma)
    if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
      warning("sigma is numerically not positive definite")
    }
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
  }
  else if (method == "chol") {
    retval <- chol(sigma, pivot = TRUE)
    o <- order(attr(retval, "pivot"))
    retval <- retval[, o]
  }
  lat <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994)
  lat[, (ncol(sigma)-ncol(latent)):ncol(sigma)] <- latent
  retval <- lat %*% retval
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}
