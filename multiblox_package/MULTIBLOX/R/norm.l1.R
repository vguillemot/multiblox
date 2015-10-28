#' Computes the l1 norm of a vector.
#' 
#' @param beta vector.
#' @return sum |beta_i|
#' @keywords internal
norm.l1 <-
function(beta){
  return(sum(abs(beta)))
}
