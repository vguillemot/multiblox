istacox.predict <- function(model.train, x.test, y.test, lambda){
  #   model.train beta estimates from the train set, 
  #   x.test matrix of covariates of the test set, 
  #   y.test outcome of the test set, 
  
  library(survival)
  source("/home/philippe/github/multiblox/istacox_MapRedR/scripts/sparse.partial.loglik.R")
  
  pll <- sparse.partial.loglik(model=model.train, newdata=x.test, newy=y.test, lambda=lambda)
  
  return(est=pll)
  
}