istacox.lambda.tune <-
  function(X, y, trainmat, i, outer_it, lambda.grid, scale=T, method="CV", metric="exactConcordanceIndex", adaptative=TRUE, fast=TRUE){
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
  source("istacox.R")
  source("istacox.predict.R")
  source("istacox.score.R")
  #source("functions.R")
    
    N <- nrow(X) # nb of individuals
    
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
    x.o <- X.train[order(y.train[, 1]), ]
    y.o <- as.data.frame(y.train[order(y.train[, 1]), ])

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
      
      ##2) Estimating the betas on the training set x, I, R, D, alpha, gamma, eps = 0.001, max.iter = 10000, beta.init = NULL
      res[[l]]<- istacox(x=x.o, I.train, R.train, alpha=0.5*lambda.grid[l,], gamma=0.25*lambda.grid[l,], beta.init=beta0, max.iter=100, adaptative=adaptative, fast=fast)
      beta.train <- res[[l]]$beta

      ### predict
      pred <- istacox.predict(model.train=res[[l]], x.train=X.train, y.train=y.train, x.test=X.test)      
      
      ##5) get and accumulate the score
      pred.score[[l]] <- istacox.score(pred$ychapo, as.matrix(y.test), metric=metric)$perf
      beta0 <- res[[l]]$beta # update beta0 for warm restart
    }
    return(list(res=res, pred.score=pred.score))
  }
