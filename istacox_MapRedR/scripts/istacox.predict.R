istacox.predict <- function(model, x, y, lambda, type=c("spll", "deviance", "pi")){
  #   model.train beta estimates from the train set, 
  #   x.test matrix of covariates of the test set, 
  #   y.test outcome of the test set, 
  #   for multiblox, this function will compute the partial loglikelihood for the training set
  
  library(survival)
  source("/home/philippe/github/multiblox/istacox_MapRedR/scripts/sparse.partial.loglik.R")
  #source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/sparse.partial.loglik.R")
  
  null_model <- rep(0, length(model))
  
  if(type=="spll"){
    res <- sparse.partial.loglik(model=model, newdata=x, newy=y, lambda=lambda)
  } else if(type=="deviance"){
    res <- -2*(sparse.partial.loglik(model=model, newdata=x, newy=y, lambda=lambda) - sparse.partial.loglik(model=null_model, newdata=x, newy=y, lambda=lambda))
  } else if(type=="pi"){
    res <- 
  } else {
    print("Sorry, this type of score is not yet implemented !")
  }
  
  
  return(est=pll)
  
}