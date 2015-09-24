istacox.lambda.tune <-
  function(X, y, trainmat, i, outer_it, lambda.grid, scale=T, method="CV", metric="spll", adaptative=TRUE, fast=TRUE){
    # x is a list of B (p_k by n) matrices of predictors
    # y dataframe of survival times and censoring, 
    # trainmat matrix of training sets indices
    # outer_it outerloop iteration
    # lambda.grid is a matrix containing all combinations of lambdas over blocks to test
    # adaptative : if true, compute an adaptative step for each ista iteration
    # TODO : propose several different metrics : inner_cv_metric=c("auc", "F", "wrec")
    
    cat("### outer CV : ", outer_it , ", inner CV iteration : ", i, "\n")
    
    library(MASS)
    library(CMA)
    library(survival)
  #  library(multiblox)
  source("/home/philippe/github/multiblox/istacox_MapRedR/scripts/istacox.R")
  source("/home/philippe/github/multiblox/istacox_MapRedR/scripts/istacox.predict.R")
  source("/home/philippe/github/multiblox/istacox_MapRedR/scripts/istacox.score.R")
  source("/home/philippe/github/multiblox/istacox_MapRedR/scripts/functions.R")
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
    
    beta.train <- eta.train <- eta.test <- NULL
    pred.score <- res <- model <- cox.model <- NULL
    
    ### Beta initialization without any intercept
    beta0 <- NULL
    
    ##1) get train and test
    ind <- trainmat[i,]
    # indexing
    y.train <- y[ind,]
    y.test <- y[-ind,]
#     print(y.test)
#     X.train <- lapply(X, function(mm) mm[ind, ])
#     if (method=="LOOCV") {
#       X.test <- lapply(X, function(mm) t(as.matrix(mm[-ind,])))
#     } else {
#       X.test <- lapply(X, function(mm) { as.matrix(mm[-ind,])})
#     }
#     # scaling
#     if (scale) {
#       X.train <- lapply(X.train, function(mm) scale2(mm))
#       scl_fun <- function(data, scaled) {
#         scale(data, center = attr(scaled, "scaled:center"),
#               scale = attr(scaled, "scaled:scale")) }
#       X.test <- mapply(scl_fun, X.test, X.train, SIMPLIFY=FALSE)
#     } 

    ### pour un bloc
    X.train <- X[[1]][ind, ]
    if (method=="LOOCV") {
      X.test <- t(as.matrix(X[[1]][-ind,]))
    } else {
      X.test <- as.matrix(X[[1]][-ind,])
    }
    # scaling
    if (scale) {
      X.train <- scale2(X.train)
      scl_fun <- function(data, scaled) {
        scale(data, center = attr(scaled, "scaled:center"),
              scale = attr(scaled, "scaled:scale")) }
      X.test <- scl_fun(X.test, X.train)
    } 


    ### Sets of patient at risk at Ti
   
    x.o <- X.train[order(y.train[, 1]), ]
    
    y.o <- as.data.frame(y.train[order(y.train[, 1]), ])

    # uncensored patients
    I.train <- which(y.o$status==1)
    # patients at risk
    R.train <- lapply( which(y.o$status==1) , function(i) which( y.o$time >= y.o$time[i] ) )
    names(R.train) <- paste0("R", which(y.o$status==1))


    ### pathwise warm restart
### en monobloc, la grille de lambda n'est plus en deux dimensions...
    if(B==1){
      lambda.grid <- t(t(lambda.grid))
    }
    for (l in 1:length(lambda.grid)){
      cat("lambda : ", l, " -> ")
#       for (bid in 1:length(lambda.grid[l,])) {
#         cat(" ", lambda.grid[l,bid], sep="")
#       }
#       cat("\n")
      
      ##2) Estimating the betas on the training set x, I, R, alpha, gamma, eps = 0.001, kmax = 10000
      print(adaptative)
      res[[l]]<- istacox(X=x.o, I.train, R.train, alpha=0.5*lambda.grid[l], gamma=0.25*lambda.grid[l], kmax=1000, ada=as.logical(adaptative), fast=as.logical(fast))
      #beta.train <- res[[l]]$beta

      ### predict
      pred <- istacox.predict(model=res[[l]], x.test=X.train, y.test=y.train, lambda=lambda.grid[l])
      
      ##5) get and accumulate the score
      ## formultiblox score is the partial loglikelihood for the training test
      ## the CV loglikelihood will be computed by the inner reducer
      pred.score[[l]] <- istacox.score(as.matrix(y.test), pred$est, type=metric)$perf
      beta0 <- res[[l]]$beta # update beta0 for warm restart
    }
    return(list(res=res, pred.score=pred.score))
  }
