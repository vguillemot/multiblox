#' Gradient
#' 
#' @param X current data block.
#' @param beta current beta.
#' @param I list of non censored individuals, ranked by event time.
#' @param R list of sets of individuals at risk at each non censored ranked time.
#' @param gamma L2-norm shrinkage parameter.
#' @param link is the total link between the current block and the other blocks.
#' @return gradient of Cox partial likelihood
#' @keywords internal
grad <- function(X, beta, I, R, gamma, link) {
  wij <- mapply( function(i, j) t(exp( X[j, ]%*%beta - log(sum( exp(X[R[[sprintf("R%i", i)]], ]%*%beta))))), I, R)
#   N <- length(I)
#   wij <- matrix(NA, N, N) 
#   for (i in 1:N) {
#     for (j in 1:N) {
#       wij[i,j] <- exp( X[j,,drop=F] %*% beta - log(sum(exp(X[R[[i]],,drop=F] %*% beta))) )
#     }
#   }
#   xj_prim_beta <- lapply(I, function(j)
#   wij_num_log <-  function(r) log(sum( exp(X[R[[sprintf("R%i", i)]], ]%*%beta))))), I, R)
#   wij <- exp( xj_prim_beta - sum(xj_prim_beta) )
  names(wij) <- names(R)
  # xbar <- t(sapply(names(R), function(r) as.matrix(wij[[r]])%*%as.matrix(X[R[[r]], ])) )
  xbar <- t(sapply(names(R), function(r) as.matrix(wij[[r]])%*%X[R[[r]],,drop=FALSE]) )
  # xbar <- t(sapply(names(R), function(r) wij[,R[[r]]]%*%X[R[[r]],]) )
  grad <- - colSums( X[I, ] - xbar ) + gamma*beta - link
}
