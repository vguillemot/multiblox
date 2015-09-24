istacox_partial_loglik <- function(x, I, R, D, beta, lambda){
  # multiblock version
#   x covariates matrix ordered by survival times, 
#   I ordered indices of non-censored patients, 
#   R sets of patients at risk for each non-censored times,
#   D design matrix, 
#   beta coefficients of the multiblox model, 
#   lambda shrinkage parameter.
  
  B <- length(x) #number of blocks
  n <- nrow(x[[1]]) #number of individuals
  p <- sapply(x, ncol) #number of covariates in each block
  
  lp <- S <- shr <- L <- NULL
   
  for (b in 1:B){
    ### loglikelihood term
    lp[[b]] <- mapply( function(i, j) x[[b]][ j, ]%*%beta[[b]] - log( sum( exp(x[[b]][R[[sprintf("R%i", i)]], ]%*%beta[[b]]))), I, R)
    names(lp[[b]]) <- names(R)
    
    ### covariance term 
    S[[b]] <- rep(0, p[[b]])
    for(d in setdiff(1:B, b)){
      S[[b]] <- S[[b]] + D[b, d] * t(beta[[b]])%*%t(x[[b]])%*%x[[d]]%*%beta[[d]]
    }
    
    ### shrinkage term
    shr[[b]] <- 0.5*lambda[[b]] * ((norm(beta[[b]], 2))^2 - 1)
    
    L[[b]] <- sum(lp[[b]]) + S[[b]] - shr[[b]]
  }
  
  LB <- lapply(L, sum)
  
  return(L = LB)
}