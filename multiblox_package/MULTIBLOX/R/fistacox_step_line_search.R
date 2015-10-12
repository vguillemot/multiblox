fistacox_step_line_search <-
function(X, u, I, R, t=10, tau=0.95, alpha, link, kmax=1000){
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
