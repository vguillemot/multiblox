istacox <-
function(X, I, R, alpha, kmax=1000, epsilon=1e-4, 
                    fast=FALSE, ada=FALSE, link, beta_init) {
  ### X is now a list of B matrices (blocks)
  ### b is the block to be treated by istacox
  ### D is a B by B matrix indicating which blocks are connected to each other
  ### I list of non censored individuals, ranked by event time.
  ### R list of sets of individuals at risk at each non censored ranked time
  ### alpha L1-norm shrinkage parameter
    
  p <- ncol(X)
  n <- nrow(X)
  
  if (!is.null(beta_init)){
    betaold <- beta_init
  } else {
    betaold <- rnorm(p) 
  }
  
  t <- 1/max(eigen(t(X)%*%X)$values)
  # print("1er grad")
  betanew <- prox(betaold - t*grad(X, betaold, I, R, alpha, link),t,alpha)
  
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
    betanew <- prox(u - t*grad(X, u, I, R, alpha, link), t, alpha)
#     print("arret ?")
#     print(k)
# print(betaold)
# print(betanew)
    if (sum((betanew-betaold)**2) < epsilon ) break
  }
  return(list(beta=betanew, k=k))
}
