istacox.predict <- function(model, x, y, D, lambda, type=c("spll", "deviance", "pi")){
  #   model list of B beta estimates from the train set, 
  #   x list of B matrices of covariates of the test set, 
  #   y list of B outcomes of the test set,  
  #   D B by B matrix indicating the links between blocks
  #   spll = sparse partial loglikelihood
  #   deviance = related to the null model
  #   pi = prognostic index (cf Bovelstad 2007)
  
  library(survival)
  source("/home/philippe/github/multiblox/istacox_MapRedR/scripts/sparse.partial.loglik.R")
  source("/home/philippe/github/multiblox/istacox_MapRedR/scripts/partial.loglik.R")
  source("/home/philippe/github/multiblox/istacox_MapRedR/scripts/link.R")
  # source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/sparse.partial.loglik.R")
  
  B <- length(x)
  
  spll <- res <- NULL
  
  null_model <- list()
  null_model[["beta"]] <- rep(0, length(model[["beta"]]))
   
  total_link <- NULL
  for (b in 1:B){
    link[b] <- link(X=x, D=D, b=b, beta.init=model)
  }
  total_link <- sum(link)
  
  if(type=="spll"){
    print("sparse loglik")
    for (b in 1:B){
      spll[[b]] <- sparse.partial.loglik(model=model[[b]], newdata=x[[b]], newy=y, lambda=lambda[b])
    }
    res <- sum(unlist(spll)) - total_link
  } else if(type=="deviance"){
    print("deviance")
    for (b in 1:B){
#       pll[[b]] <- sparse.partial.loglik(model=model[[b]], newdata=x[[b]], newy=y, lambda=lambda[b])
#       dev_null[[b]] <- sparse.partial.loglik(model=null_model, newdata=x[[b]], newy=y, lambda=lambda[b])
      pll[[b]] <- partial.loglik(model=model[[b]], newdata=x[[b]], newy=y)
      dev_null[[b]] <- partial.loglik(model=null_model, newdata=x[[b]], newy=y)
    }
    total_pll <- sum(unlist(pll)) - total_link
    print(paste("Loglik du modele : ", total_pll, sep=""))
    total_null <- sum(unlist(dev_null))
    print(paste("Loglik du modele nul : ", total_null, sep=""))
    res <- -2*(total_spll - total_null)
    print(paste("Deviance du modele par rapport au modele mul : ", res, sep=""))
  } else if(type=="pi"){
    print("pronostic index")
    for (b in 1:B){
      dat[[b]] <- x[[b]]%*%model[[b]]
    }
    lp <- Reduce(cbind, x = dat)
    fit <- lapply(coxph(Surv(y[, 1], y[,2])~lp)
    res <- summary(fit)$logtest["pvalue"]
  } else {
    print("Sorry, this type of score is not yet implemented !")
  }

  return(list(est=res))
  
}