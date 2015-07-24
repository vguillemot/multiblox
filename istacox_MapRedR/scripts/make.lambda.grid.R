make.lambda.grid <-
  function(X, path=c("naive", "smart", "norm1")){
    ##################################################################
    ###
    ###  Ridge Lambda grid definition cf Coxnet, Simon et al, JSS 2011.
    ###
    ##################################################################

    B <- length(X) # nb of blocks
    
    # scaling
    scale <- TRUE
    if (scale) {
      X <- lapply(X, function(mm) scale2(mm))
    }
    
    lambda.max <- rep(0, B)
    eps <- 0.005
    for (i in 1:B){
      lambda.max[[i]] <- 0.5*sum(apply(X[[i]], 2, function(c) norm(c, "2")))
    }
    lambda.min <- eps*lambda.max
    
    # -1 because sequence starts at 0 (0:nb.lambda) and -1 to include huge lambda (<Inf)
    nb.lambda <- 10 - 1 - 1 # corresponds to m, in the paper
    l <- NULL
    # 0 to generate keep space for huge lambda
    j <- c(0, 0:nb.lambda) 
    for(i in 1:B){
      l[[i]]<- lambda.max[[i]]*(lambda.min[[i]]/lambda.max[[i]])^(j/nb.lambda)
      # introduce huge lambda at the beginning
      l[[i]][[1]] <- 1e+20
    }
    if(B==1){
      print("un seul bloc !!")
      lg <- as.matrix(expand.grid(l))
      lambda.grid <- matrix(0, nrow=length(lg), ncol=1)
      lambda.grid[, 1] <- lg
    }else{
      lambda.grid <- as.matrix(expand.grid(l))
    }
    
    
    ### pathwise : naive ranking
    if (path == "naive") {
      n1 <- nrow(lambda.grid):1
    } else if (path  == "norm1") {
      n1 = apply(lambda.grid, 1, function(u) u[[1]] + u[[2]]) #seulement pour B=2
    } else if (path == "smart") {
      ## TODO ##
      n1 <- nrow(lambda.grid):1
    }
    #lambda.grid = lambda.grid[order(n1, decreasing=T),]
    return(lambda.grid[order(n1, decreasing=T),])
  }

make.grid <-
  function(data){
  return (make.lambda.grid(data$X, path="naive"))   
}