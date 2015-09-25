istacox.predict <- function(model, x, y, lambda, type=c("spll", "deviance", "pi")){
  #   model.train beta estimates from the train set, 
  #   x.test matrix of covariates of the test set, 
  #   y.test outcome of the test set,  
  #   spll = sparse partial loglikelihood
  #   deviance = related to the null model
  #   pi = prognostic index (cf Bovelstad 2007)
  
  library(survival)
  source("/home/philippe/github/multiblox/istacox_MapRedR/scripts/sparse.partial.loglik.R")
  #source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/sparse.partial.loglik.R")
  
  null_model <- rep(0, length(model))
  
  if(type=="spll"){
    res <- sparse.partial.loglik(model=model, newdata=x, newy=y, lambda=lambda)
  } else if(type=="deviance"){
    res <- -2*(sparse.partial.loglik(model=model, newdata=x, newy=y, lambda=lambda) - sparse.partial.loglik(model=null_model, newdata=x, newy=y, lambda=lambda))
  } else if(type=="pi"){
    dat <- list(pi=t(x)%*%model$beta, time=y[, 1], status=y[,2])
    fit <- coxph(Surv(time, status)~pi, data=dat)
    res <- summary(fit)$logtest["pvalue"]
  } else {
    print("Sorry, this type of score is not yet implemented !")
  }

  return(est=res)
  
}