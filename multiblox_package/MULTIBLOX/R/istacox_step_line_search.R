istacox_step_line_search <-
function(X, beta, I, R, t=10, tau=0.5, alpha, link, kmax=1000){
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
