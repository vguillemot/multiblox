#' Compute the link between blocks.
#' 
#' @param X a list of B matrices
#' @param D is a B by B binary design matrix describing connections between blocks
#' @param b is the block to be treated by istacox
#' @param beta list of B cox coefficient vectors
#' @return The sum of covariance between linked blocks
link <- function(X, D, b, beta.init){
  
  B <- length(X)
  p <- lapply(X, ncol)
  
  S <- rep(0, p[[b]])
  
  if (!is.null(beta.init)){
    beta <- beta.init
  }else{
    f <- function(z) z/1e2
    for (c in 1:B){
      beta[[c]] <- rep(0, p[[c]])
    }
  }
  eta <- mapply(X, beta, FUN=function(a, b) a%*%b, SIMPLIFY=FALSE)  
  for(d in setdiff(1:B, b)){
    S <- S + D[b, d] * t(X[[b]])%*%eta[[d]]
  }
  return(S)
}
