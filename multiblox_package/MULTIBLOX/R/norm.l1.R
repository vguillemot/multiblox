#' Computes the L1-norm of a vector
#' 
#' @param beta is a vector of values
#' @return the L1-norm of beta = sum(|beta|)
#' @keywords internal
norm.l1 <-
function(beta){
  return(sum(abs(beta)))
}
