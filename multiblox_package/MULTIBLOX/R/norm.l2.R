#' Computes the L2-norm of a vector
#' 
#' @param beta is a vector of values
#' @return the L2-norm of beta = square root (sum(|beta|))
#' @keywords internal
norm.l2 <-
function(beta){
  return(sqrt(sum(beta^2)))
}
