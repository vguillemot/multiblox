generate_monoblock_survival_data <-
function(n, p, p_frac){
  set.seed(10101)
  N <- n
  nzc <- p*p_frac
  x <- matrix(rnorm(N*p), nrow=N, ncol=p)
  betastar <- rnorm(nzc)
  fx <- x[,seq(nzc)]%*%(betastar*p_frac)
  h0 <- 1
  hx <- h0*exp(fx)
  ty <- rexp(N,hx)
  # censoring indicator
  tcens <- rbinom(n=N, prob=.3, size=1)
  # y=Surv(ty,1-tcens) with library(survival)
  y <- cbind(time=ty,status=1-tcens)
 
  X=list()
  X[[1]] = x
  my_var_names <- sprintf("blk1var%d", 1:p)
  colnames(X[[1]]) <- my_var_names
  
  return (list(X=X, y=as.matrix(y), my_var_names=my_var_names, h0=h0))
  
}
