#' Computes the negative step that updates the proximal gradient in adaptative ISTA
#' 
#' @param X is a matrix of observations
#' @param t is the update step
#' @param beta are the current coefficients estimates
#' @param alpha is the sparsity parameter
#' @param g is the gradient
#' @return the negative step to update the proximal gradient
#' @keywords internal
neg_step <-
function(X, t, beta, alpha, g){
  return((1/t) * (beta - prox(beta - t * g, t, alpha)))
}
