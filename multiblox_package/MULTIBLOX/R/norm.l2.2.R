#' Computes the squared L2-norm of a vector
#' 
#' @param beta is a vector of values
#' @return the squared L2-norm of beta = sum(|beta|^2)
#' @keywords internal
norm.l2.2 <-
function(beta){
  return(sum(beta^2))
}
