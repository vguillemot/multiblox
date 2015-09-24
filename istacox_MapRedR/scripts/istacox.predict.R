istacox.predict <- function(model.train, x.train, y.train, lambda, type=c("spll", "deviance", "logrank")){
  #   model.train beta estimates from the train set, 
  #   x.test matrix of covariates of the test set, 
  #   y.test outcome of the test set, 
  #   for multiblox, this function will compute the partial loglikelihood for the training set
  
  library(survival)
  source("/home/philippe/github/multiblox/istacox_MapRedR/scripts/sparse.partial.loglik.R")
  #source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/sparse.partial.loglik.R")
  
  if(type=="spll"){
    res <- sparse.partial.loglik(model=model.train, newdata=x.train, newy=y.train, lambda=lambda)
  } else if(type=="deviance"){
    res <- 
  } else if(type=="logrank"){
    res <- 
  } else {
    print("Sorry, this type of score is not yet implemented !")
  }
  
  
  return(est=pll)
  
}