
# grad <- function(X, beta, I, R, gamma) {
#   wij <- mapply( function(i, j) t(exp( X[j, ]%*%beta) / sum( exp(X[R[[sprintf("R%i", i)]], ]%*%beta))), I, R)
#   names(wij) <- names(R)
#   xbar <- t(sapply(names(R), function(r) wij[[r]]%*%X[R[[r]], ]) )
#   grad <- -colSums( X[I, ] - xbar ) + gamma*beta 
# }


grad <- function(X, beta, I, R, alpha, link) {
#   print(X)
#   print(beta)
#   print(X[1,]%*%beta)
#   print(I)
  # print(R)
  # wij <- mapply( function(i, j) t(exp(matrix(X[j, ], nrow=1)%*%matrix(beta, ncol=1)) / sum( exp(X[R[[sprintf("R%i", i)]], ]%*%beta))), I, R, SIMPLIFY = FALSE)
  # wij <- mapply( function(i, j) exp(beta%*%X[j,] / sum( exp(beta%*%X[R[[sprintf("R%i", i)]], ]))), I, R, SIMPLIFY = FALSE)
  wij <- mapply( function(i, j) t(exp( X[j, ]%*%beta) / sum( exp(X[R[[sprintf("R%i", i)]], ]%*%beta))), I, R)
#   print(dim(wij[[1]]))
#   print(dim(as.matrix(X[R[[1]], ])))
  names(wij) <- names(R)
  xbar <- t(sapply(names(R), function(r) as.matrix(wij[[r]])%*%as.matrix(X[R[[r]], ])) )
  grad <- - colSums( X[I, ] - xbar ) + (1-alpha)*beta - link
}

neg_step <- function(X, t, beta, alpha, g){
  return((1/t) * (beta - prox(beta - t * g, t, alpha)))
}

R_ell <- function(X, beta, I, R, alpha){
  truc <- mapply( function(i, j) {matrix(X[j, ], nrow=1)%*%matrix(beta, ncol=1) -log(sum( exp(X[R[[sprintf("R%i", i)]], ]%*%beta)))}, I, R, SIMPLIFY = F)
  - sum(unlist(truc)) + (1-alpha)/2 * norm.l2.2(beta)
}

prox <- function(beta,t,alpha) ifelse(abs(beta)<t*alpha,0,beta-t*alpha*sign(beta))

istacox_step_line_search <- function(X, beta, I, R, t=10, tau=0.5, alpha, link, kmax=1000){
  R_l <- R_ell(X, beta, I, R, alpha)
  g <- grad(X, beta, I, R, alpha, link)
  Gt <- neg_step(X, t, beta, alpha, g)
  Rt_l <- R_ell(X, beta - t*Gt, I, R, alpha)

  for (k in 1:kmax) {
    if (Rt_l <= R_l -t*t(g) %*% Gt + t/2*sum(Gt**2)) break
    t <- tau*t
    Gt <- neg_step(X, t, beta, alpha, g)
    Rt_l <- R_ell(X, beta - t*Gt, I, R, alpha)
  }
  return(t)
}


fistacox_step_line_search <- function(X, u, I, R, t=10, tau=0.95, alpha, link, kmax=1000){
  Ru_l <- R_ell(X, u, I, R, alpha)
  gradu <- grad(X, u, I, R, alpha, link)
  x <- prox(u - t*gradu, t, alpha)
  Rx_l <- R_ell(X, x, I, R, alpha)

  for (k in 1:kmax) {
    if (Rx_l <= Ru_l +  crossprod(gradu, x-u) + 1/(2*t)*sum((x-u)**2)) break
    t <- tau*t
    x <- prox(u - t*gradu,t, alpha)
    Rx_l <- R_ell(X, x, I, R, alpha)
  }
  return(t)
}

istacox <- function(X, I, R, alpha, kmax=1000, epsilon=1e-4, 
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
