link <- function(X, D, b, beta.init){
  ### X a list of B matrices
  ### D is a B by B binary design matrix describing connections between blocks
  ### b is the block to be treated by istacox
  ### beta list of B cox coefficient vectors
  
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


relax_multiblox <- function(x, I, R, D, lambda = 0, eps = 0.001, max.iter = 10000, beta.init = NULL){
  ### Block relaxation for multiblox
  
    # x is a list of B (p_k by n) matrices of predictors
    # I vector of indexes of uncensored patients, 
    # R list of vector of indexes of patients at risk for each ordered time Ti, 
    # D is a B by B binary design matrix describing connections between blocks
    # lambda.opt is a vector of shrinkage parameters, one for each block
    # eps is a real corresponding to the tolerance for convergence
    # max.iter is an integer corresponding to the maximum nb of ALS iterations
    # beta.init is a p by 1 vector of the initial values of beta
    
    source("istacox.R")
    
    B <- length(x) #number of blocks
    n <- nrow(x[[1]])
    p <- sapply(x, ncol) #number of covariates in each block
    
    beta <- beta_new <- e <- S <- eta <- cvg <- NULL
    
    ### Initialization step
    iter <- 1
    cvg[[iter]] <- 10
    beta <- list()
    
    if (!is.null(beta.init)){
      beta <- beta_new <- beta.init
    }else{
      f <- function(x) x/1e2
      for (c in 1:B){
        beta[[c]] <- rep(0, p[c])
      }
      beta_new <- beta
    }
    iter.inner <- 0
    div.inner <- 0
    max.iter.inner <- 200
    consec_max.iter <- 0
    beta_new <- beta
    
    for (iter in 1:max.iter){
      for (b in 1:B) {
        
        link <- link(X, D, b, beta)
        istacox_res <- istacox(x=x, I=I, R=R, D=D, b=b, alpha=0.5*lambda, kmax=1000, epsilon=1e-6, 
                           fast=FALSE, ada=FALSE, link=link, beta_init=beta)
        iter.inner <- iter.inner + istacox_res$k
        print(iter)
        print(iter.inner)
        if (max.iter.inner == istacox_res$k) {
          div.inner <- div.inner + 1
          consec_max.iter <- consec_max.iter + 1
        } else {
          consec_max.iter <- 0
        }
        #       beta_new <- lapply(NR_step$beta, beta_norm2)
        beta_new <- istacox_res$beta
      }
      d <- mapply("-", beta, beta_new, SIMPLIFY=FALSE)
      e <- sapply(d, base::norm, "f")
      if(max(e)<eps) break
      beta <- beta_new
      if (consec_max.iter >= (2 * B)) break # all the blocks diverge 2 times
    }
    #eta <- mapply(x, beta, FUN=function(a, b) a%*%b, SIMPLIFY=FALSE) ### bug 3 blocs
    print(paste("Block relaxation iter (inner - div): ", iter , "(", iter.inner, "-", div.inner, ")" ))
    return(list(beta = beta_new, convergence = unlist(e), niter=iter))
  }
  
}