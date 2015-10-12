grad <-
function(X, beta, I, R, alpha, link) {
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
