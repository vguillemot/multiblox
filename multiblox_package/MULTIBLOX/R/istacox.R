#' Computes the beta coefficients for a sparse cox model via ISTA
#' 
#' @param X is a list of B matrices (blocks)
#' @param I list of non censored individuals, ranked by event time.
#' @param R list of sets of individuals at risk at each non censored ranked time
#' @param alpha L1-norm shrinkage parameter
#' @param kmax is the maximal number of iterations.
#' @param epsilon is the required precision.
#' @param fast is a boolean used to specify if FISTA (TRUE) is used instead of ISTA (FALSE, default).
#' @param ada is a boolean used to specify is the step must be chosen at each iteration (TRUE) or not (FALSE, default).
#' @param link is the total link between the current block and the other blocks.
#' @param beta_init is the initial value of the weight vector (used mainly for warm restarts).
#' @return beta coefficients that maximize the Cox sparse partial likekihood, and k number of iterations
istacox <- function(X, I, R, alpha, kmax=1000, epsilon=1e-4, 
                    fast=FALSE, ada=FALSE, link, beta_init) {

    
  p <- ncol(X)
  n <- nrow(X)
  
  if (!is.null(beta_init)){
    betaold <- beta_init
  } else {
    betaold <- rnorm(p) 
  }
  
  t <- 1/max(eigen(t(X)%*%X)$values)
  # print("1er grad")
  betanew <- prox(betaold - t*grad(X, betaold, I, R, alpha, link),t,alpha)
  
  for (k in 2:kmax) {
    if (fast) {
      u <- betanew + (k-1) / (k+2) * (betanew - betaold)
    } else {
      u <- betanew
    }
    if (ada) {
      if (fast) {
        t <- fistacox_step_line_search(X, u, I, R, t, tau=0.95, alpha, link)
      } else {
        t <- istacox_step_line_search(X, betaold, I, R, t, tau=0.95, alpha, link)
      }
    }
    betaold <- betanew
    betanew <- prox(u - t*grad(X, u, I, R, alpha, link), t, alpha)
    if (sum((betanew-betaold)**2) < epsilon ) break
  }
  return(list(beta=betanew, k=k))
}
