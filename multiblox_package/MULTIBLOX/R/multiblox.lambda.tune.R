#' Computes all sparse cox models for a grid of parameter lambda
#' 
#' @param X is a list of B matrices (blocks)
#' @param y is the survival variable including time and status
#' @param D is the matrix of design storing the links among blocks
#' @param trainmat a matrix of observations indices indicating the training set for each CV fold of model selection.
#' @param i iteration of inner CV loop
#' @param outer_it iteration of outer CV loop
#' @param lambda.grid a grid of tuples of lambda parameters
#' @param scale a boolean indicating if scaling should be performed
#' @param method is either "CV", "LOOCV" or "MCCV" to choose CV method
#' @param metric is either "spll" for sparse partial loglikelihood, "deviance" or "pi" for Prognostic index is computed for model assessment and selection.
#' @param adaptative boolean used to specify is the step must be chosen at each iteration (TRUE) or not (FALSE, default).
#' @param fast boolean used to specify if FISTA (TRUE) is used instead of ISTA (FALSE, default).
#' @param path is either "naive", "smart" (not implemented yet) or "norm1" (only if B=2)
#' @return lambda.grid, 10 lambda values per block ie 10^B tuples of lambdas
multiblox.lambda.tune <-
function(X, y, D, trainmat, i, outer_it, lambda.grid, scale=T, method="CV", metric="spll", adaptative=TRUE, fast=TRUE){
    # x is a list of B (p_k by n) matrices of predictors
    # y dataframe of survival times and censoring, 
    # D is a B by B binary design matrix describing connections between blocks
    # trainmat matrix of training sets indices
    # outer_it outerloop iteration
    # lambda.grid is a matrix containing all combinations of lambdas over blocks to test
    # adaptative : if true, compute an adaptative step for each ista iteration
    # TODO : propose several different metrics : inner_cv_metric=c("auc", "F", "wrec")
    
    cat("### outer CV : ", outer_it , ", inner CV iteration : ", i, "\n")
    
    library(MASS)
    library(CMA)
    library(survival)
    library(MULTIBLOX)
#   source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/relax_multiblox.R")
#   source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/istacox.predict.R")
#   source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/istacox.score.R")
#   source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/functions.R")
#   source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/istacox.R")
#   source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/istacox.predict.R")
#   source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/istacox.score.R")
#   source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/functions.R")
    
    
    B <- length(X) # nb of blocks
  if(B==1){
    N <- nrow(X[[1]]) # nb of individuals
  }else{
    N <- nrow(X) # nb of individuals
  }
print(D)
      
      beta.train <- eta.train <- eta.test <- NULL
      pred.score <- res <- model <- cox.model <- cv <- NULL
      
      ### Beta initialization without any intercept
      beta0 <- NULL
      for (b in 1:B){
        beta0[[b]] <- matrix(0, nrow=(ncol(X[[b]])), ncol=1)
      }
      
      ##1) get train and test
      ind <- trainmat[i,]
      # indexing
      y.train <- y[ind,]
      y.test <- y[-ind,]
      #     print(y.test)
      X.train <- lapply(X, function(mm) mm[ind, ])
      if (method=="LOOCV") {
        X.test <- lapply(X, function(mm) t(as.matrix(mm[-ind,])))
      } else {
        X.test <- lapply(X, function(mm) { as.matrix(mm[-ind,])})
      }
      # scaling
      if (scale) {
        X.train <- lapply(X.train, function(mm) scale2(mm))
        scl_fun <- function(data, scaled) {
          scale(data, center = attr(scaled, "scaled:center"),
                scale = attr(scaled, "scaled:scale")) }
        X.test <- mapply(scl_fun, X.test, X.train, SIMPLIFY=FALSE)
      } 
      
      ### Sets of patient at risk at Ti
      x.o <- list()
      x.o <- lapply(X.train, function(b) b[order(y.train[, 1]), ])
      y.o <- as.data.frame(y.train[order(y.train[, 1]), ])
      colnames(y.o) <- c("time", "status")
#       print(length(x.o))
#       print(dim(x.o[[1]]))
#       print(dim(x.o[[2]]))
#       print(dim(x.o[[3]]))
      
      # uncensored patients
      I.train <- which(y.o$status==1)
      # patients at risk
      R.train <- lapply( which(y.o$status==1) , function(i) which( y.o$time >= y.o$time[i] ) )
      names(R.train) <- paste0("R", which(y.o$status==1))
      
      
      ### pathwise warm restart              
      for (l in 1:nrow(lambda.grid)){
        cat("lambdas : ", l, " -> ")
        for (bid in 1:length(lambda.grid[l,])) {
          cat(" ", lambda.grid[l,bid], sep="")
        }
        cat("\n")
  
      ##2) Estimating the betas on the training set x, I, R, alpha, gamma, eps = 0.001, kmax = 10000
      print(adaptative)
      res[[l]]<- relax_multiblox(x=x.o, I.train, R.train, D=D, lambda=lambda.grid[l,], max.iter=1000, ada=as.logical(adaptative), fast=as.logical(fast), beta.init=beta0)
      beta.train[[l]] <- res[[l]]$beta

      ### predict
      pred <- istacox.predict(model=beta.train[[l]], D=D, x=x.o, y=y.o, lambda=lambda.grid[l,], type=metric)
      
      ##5) get and accumulate the score
      ## formultiblox score is the partial loglikelihood for the training test
      ## the CV loglikelihood will be computed by the inner reducer
      pred.score[[l]] <- istacox.score(as.matrix(y.test), pred[["est"]])[["perf"]]
      beta0 <- res[[l]][["beta"]] # update beta0 for warm restart
    }
    return(list(res=res, pred.score=pred.score))
  }
