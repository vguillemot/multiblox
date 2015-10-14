#' Computes the step for an iteration of ISTA Cox regression with a line search algorithm.
#' 
#' @param X current data block.
#' @param beta current beta.
#' @param I list of non censored individuals, ranked by event time.
#' @param R list of sets of individuals at risk at each non censored ranked time.
#' @param t is the maximum step for the line search.
#' @param tau is the factor for the line search.
#' @param alpha L1-norm shrinkage parameter.
#' @param link is the total link between the current block and the other blocks.
#' @param kmax is the maximal number of iterations.
#' @return step
istacox_step_line_search <- function(X, beta, I, R, t=10, tau=0.5, alpha, link, kmax=1000){
  R_l <- R_ell(X, beta, I, R, alpha)
  g <- grad(X, beta, I, R, alpha, link)
  Gt <- neg_step(X, t, beta, alpha, g)
  Rt_l <- R_ell(X, beta - t*Gt, I, R, alpha)

  for (k in 1:kmax) {
    if (Rt_l <= R_l -t*t(g) %*% Gt + t/2*sum(Gt**2)) break
    t <- tau*t
    Gt <- neg_step(X, t, beta, alpha, g)
    Rt_l <- R_ell(X, beta - t*Gt, I, R, alpha)
  }
  return(t)
}
