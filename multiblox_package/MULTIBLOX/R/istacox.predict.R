#' Predicts (sparse log-likelihood, prognostic index or deviance) from a model and test data.
#' 
#' @param model list of B beta estimates from the train set, 
#' @param x list of B matrices of covariates of the test set, 
#' @param y list of B outcomes of the test set,  
#' @param D B by B matrix indicating the links between blocks
#' @param spll = sparse partial loglikelihood
#' @param deviance = related to the null model
#' @param pi = prognostic index (cf Bovelstad 2007)
#' @return the quality of the prediction specified by the user.
istacox.predict <-
function(model, x, y, D=NULL, lambda, type=c("spll", "deviance", "pi")){

#   library(survival)
#   source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/sparse.partial.loglik.R")
#   source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/partial.loglik.R")
#   source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/link.R")
#   source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/sparse.partial.loglik.R")
  
  B <- length(x)
  
  spll <- res <- lk <- pll <- dev_null <- NULL
  
  null_model <- matrix(0, nrow=length(model), ncol=1)
  
  total_link <- 0
  if (!is.null(D)){
    for (b in 1:B){
      lk[[b]] <- link(X=x, D=D, b=b, beta.init=model)
    }
    total_link <- sum(unlist(lk))
  }
    
  if(type=="spll"){
    print("sparse loglik")
    if (!is.null(D)){
      for (b in 1:B){
        spll[[b]] <- sparse.partial.loglik(model=model[[b]], newdata=x[[b]], newy=y, lambda=lambda[b])
      }
      res <- sum(unlist(spll)) - total_link
    } else {
      spll <- sparse.partial.loglik(model=model, newdata=x, newy=y, lambda=lambda)
      res <- spll
    }
  } else if(type=="deviance"){
    print("deviance")
    if (!is.null(D)){
      # print("multiblox")
      for (b in 1:B){
#         print(model[[b]])
#         print(x[[b]])
        pll[[b]] <- partial.loglik(model=model[[b]], newdata=x[[b]], newy=y)
        dev_null[[b]] <- partial.loglik(model=matrix(0, nrow = length(model[[b]]), ncol=1), newdata=x[[b]], newy=y)
      }
      total_pll <- sum(unlist(pll)) - total_link
      # print(paste("Loglik du modele : ", total_pll, sep=""))
      total_null <- sum(unlist(dev_null))
      # print(paste("Loglik du modele nul : ", total_null, sep=""))
      res <- -2*(total_pll - total_null)
      # print(paste("Deviance du modele par rapport au modele mul : ", res, sep=""))
    } else {
      # print("coxnet")
      pll <- partial.loglik(model=model, newdata=x, newy=y)
      # print(paste("pll : ", pll, sep=""))
      dev_null <- partial.loglik(model=matrix(0, nrow = length(model), ncol=1), newdata=x, newy=y)
      # print(paste("pll_null : ", dev_null, sep=""))
      res <- -2*(pll - dev_null)
      # print(res)
    }
  } else if(type=="pi"){
    print("pronostic index")
    if (!is.null(D)){
      for (b in 1:B){
        dat[[b]] <- x[[b]]%*%model[[b]]
      }
      lp <- Reduce(cbind, x = dat)
      fit <- lapply(coxph(Surv(y[, 1], y[,2])~lp))
      res <- summary(fit)$logtest["pvalue"]
    } else {
      dat <- as.matrix(x)%*%as.matrix(model)
      fit <- coxph(Surv(y[, 1], y[,2])~dat)
      res <- summary(fit)$logtest["pvalue"]
    }
  } else {
    print("Sorry, this type of score is not yet implemented !")
  }
  return(list(est=res))
}
