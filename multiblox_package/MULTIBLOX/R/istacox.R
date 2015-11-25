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
  
  t <- 1/(max(eigen(t(X)%*%X)$values) + alpha/2)
#   told <- 1/max(eigen(t(X)%*%X)$values)
#   print(paste("Old step : ", told))
#   DD <- apply(X, 2, function(v) sum(sapply(names(R), function(r) 1/(4*n)*diff(range(v[R[[r]]]))^2 )) )
#   # t <- 1/(max(DD) + alpha/2) 
#   t <- 1/diag(DD + alpha/2) 
#   print(paste("New step : ", min(t)))
  # print("1er grad")
  betanew <- prox(betaold - t*grad(X, betaold, I, R, 0.5*alpha, link), t, alpha)
  pll <- sum(mapply( function(i) {X[i, ]%*%betanew - log(sum( exp(X[R[[sprintf("R%i", i)]],,drop=F]%*%betanew) ))}, I))
  
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
    betanew <- prox(u - t*grad(X, u, I, R, 0.5*alpha, link), t, alpha)
    if (sum((betanew-betaold)**2) < epsilon ) {
      pll <- sum(mapply( function(i) {X[i, ]%*%betanew - log(sum( exp(X[R[[sprintf("R%i", i)]],,drop=F]%*%betanew) ))}, I))
      pen <- alpha*(norm.l1(betanew) + 0.5*norm.l2.2(betanew))
      print(paste("pll : ", pll))
      print(paste("pen : ", pen))
      print(paste("step : ", t))
      break
    }
    pllnew <- sum(mapply( function(i) {X[i, ]%*%betanew - log(sum( exp(X[R[[sprintf("R%i", i)]],,drop=F]%*%betanew) ))}, I))
    if (pllnew < pll) {
      # print("Something is wrong !")
      # t <- t/10
      # print(paste("New step : ", t))
      # link <- 0.5*link
      # print(paste("New link : ", link))
      
    }
    pll <- pllnew
    pen <- alpha*(norm.l1(betanew) + 0.5*norm.l2.2(betanew))
    # print(paste("step : ", t))
    # print(paste("pll : ", pll))
    # print(paste("pen : ", pen))
    # print(paste("critere total : ", -pll + pen - link %*% betanew, sep=""))
    
  }
  return(list(beta=betanew, k=k))
}
