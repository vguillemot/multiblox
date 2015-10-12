R_ell <-
function(X, beta, I, R, alpha){
  truc <- mapply( function(i, j) {matrix(X[j, ], nrow=1)%*%matrix(beta, ncol=1) -log(sum( exp(X[R[[sprintf("R%i", i)]], ]%*%beta)))}, I, R, SIMPLIFY = F)
  - sum(unlist(truc)) + (1-alpha)/2 * norm.l2.2(beta)
}
