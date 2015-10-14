#' Gradient
#' 
#' @param X current data block.
#' @param beta current beta.
#' @param I list of non censored individuals, ranked by event time.
#' @param R list of sets of individuals at risk at each non censored ranked time.
#' @param alpha L1-norm shrinkage parameter.
#' @param link is the total link between the current block and the other blocks.
#' @return gradient
#' @keywords internal
grad <- function(X, beta, I, R, alpha, link) {
  wij <- mapply( function(i, j) t(exp( X[j, ]%*%beta) / sum( exp(X[R[[sprintf("R%i", i)]], ]%*%beta))), I, R)
  names(wij) <- names(R)
  xbar <- t(sapply(names(R), function(r) as.matrix(wij[[r]])%*%as.matrix(X[R[[r]], ])) )
  grad <- - colSums( X[I, ] - xbar ) + (alpha/2)*beta - link
}
