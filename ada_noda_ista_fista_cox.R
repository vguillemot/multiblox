
grad <- function(X, y, beta, I, R, gamma) {
  wij <- mapply( function(i, j) t(exp( X[j, ]%*%beta) / sum( exp(X[R[[sprintf("R%i", i)]], ]%*%beta))), I, R)
  names(wij) <- names(R)
  xbar <- t(sapply(names(R), function(r) wij[[r]]%*%X[R[[r]], ]) )
  grad <- -colSums( X[I, ] - xbar ) + gamma*beta 
}

prox <- function(beta,t,alpha) ifelse(abs(beta)<t*alpha,0,beta-t*alpha*sign(beta))

istacox_step_line_search <- function(X, y, beta, I, R, t=10, tau=0.5, alpha, gamma, kmax=1000){
  R_l <- R_ell(X, y, beta, I, R, gamma)
  g <- grad(X, y, beta, I, R, gamma)
  Gt <- neg_step(X, y, t, beta, alpha, g)
  Rt_l <- R_ell(X, y, beta - t*Gt, I, R, gamma)

  for (k in 1:kmax) {
    if (Rt_l <= R_l -t*t(g) %*% Gt + t/2*sum(Gt**2)) break
    t <- tau*t
    Gt <- neg_step(X, y, t, beta, alpha, g)
    Rt_l <- R_ell(X, y, beta - t*Gt, I, R, gamma)
  }
  return(t)
}


fistacox_step_line_search <- function(X, y, u, I, R, t=10, tau=0.95, alpha, gamma, kmax=1000){
  Ru_l <- R_ell(X, y, u, I, R, gamma)
  gradu <- grad(X, y, u, I, R, gamma)
  x <- prox(u - t*gradu, t, alpha)
  Rx_l <- R_ell(X, y, x, I, R, gamma)

  for (k in 1:kmax) {
    if (Rx_l <= Ru_l +  crossprod(gradu, x-u) + 1/(2*t)*sum((x-u)**2)) break
    t <- tau*t
    x <- prox(u - t*gradu,t,alpha)
    Rx_l <- R_ell(X, y, x, I, R, gamma)
  }
  return(t)
}

istacox <- function(X, y, I, R, alpha, gamma, kmax=1000, epsilon=1e-10, 
                     fast=FALSE, ada=FALSE) {
  p <- ncol(X)
  n <- nrow(X)
  betaold <- rnorm(p)
  if (!ada) t <- 1/max(eigen(t(X)%*%X)$values)

  betanew <- prox(betaold - t*grad(X, y, betaold, I, R, gamma),t,alpha)
  
  for (k in 2:kmax) {
    if (fast) {
      u <- u + (k-1) / (k+2) * (betanew - betaold)
    } else {
      u <- betaold
    }
    if (ada) {
      if (fast) {
        t <- fistacox_step_line_search(X, y, u, I, R, t, tau=0.95, alpha, gamma)
      } else {
        t <- istacox_step_line_search(X, y, betaold, I, R, t, tau=0.95, alpha, gamma)
      }
    }
    betanew <- prox(u - t*grad(X, y, u, I, R, gamma), t, alpha)
    if (sum((betanew-betaold)**2) < epsilon ) break
    betaold <- betanew
  }
  return(list(beta=betanew, k=k))
}

