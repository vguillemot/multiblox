#' Proximal operator of the L1-norm
#' 
#' @param beta is a vector of coefficients
#' @param t is the step for each iteratio of ISTA
#' @param alpha is L1-norm sparsity parameter
#' @return prox of L1-norm
#' @keywords internal
prox <-
function(beta,t,alpha) ifelse(abs(beta)<t*alpha,0,beta-t*alpha*sign(beta))
