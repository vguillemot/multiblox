neg_step <-
function(X, t, beta, alpha, g){
  return((1/t) * (beta - prox(beta - t * g, t, alpha)))
}
