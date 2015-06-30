istacox.predict <- function(model.train, x.test, y.test){
  #   model.train beta estimates from the train set, 
  #   x.test matrix of covariates of the test set, 
  #   y.test outcome of the test set, 
  
  library(survival)
  pll <- partial_loglik(model=model.train, newdata=x.test, newy=y.test)
  
  return(est=pll)
  
}