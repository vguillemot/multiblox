#' Computes the residuals of shrinked cox regression
#' 
#' @param X is a list of B matrices (blocks)
#' @param beta is vector of current cox coefficients
#' @param I list of non censored individuals, ranked by event time.
#' @param R list of sets of individuals at risk at each non censored ranked time
#' @param alpha is the shrinkage parameter
#' @return the residuals of sparse cox regression
#' @keywords internal
R_ell <- function(X, beta, I, R, alpha){
  
  truc <- mapply( function(i, j) {X[j,,drop=FALSE] %*% beta -log(sum( exp(X[R[[sprintf("R%i", i)]],,drop=F]%*%beta)))}, I, R, SIMPLIFY = F)
  - sum(unlist(truc)) + (1-alpha)/2 * norm.l2.2(beta)
}
