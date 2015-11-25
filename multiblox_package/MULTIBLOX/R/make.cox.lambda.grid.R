#' Computes a grid of parameter for Cox model
#' 
#' @param X is a list of B matrices (blocks)
#' @param y is the survival variable, made of time and status
#' @param nl is the number of lambdas for each block in the grid
#' @param path is either "naive", "smart" (not implemented yet) or "norm1" (only if B=2)
#' @return lambda.grid, 10 lambda values per block ie 10^B tuples of lambdas
make.cox.lambda.grid <-
  function(X, y, nl=10, path=c("naive", "smart", "norm1")){
    ##################################################################
    ###
    ###  ElasticNet Lambda grid definition cf Yang et al (2013) the cocktail algorithm.
    ###
    ##################################################################
    
    B <- length(X) # nb of blocks
    n <- dim(X[[1]])[1]
    alpha <- 0.5
    
    # scaling
    scale <- TRUE
    if (scale) {
      X <- lapply(X, function(mm) scale2(mm))
    }
    
    # Ordering data
    X.o <- list()
    X.o <- lapply(X, function(b) b[order(y[, 1]), ])
    y.o <- as.data.frame(y[order(y[, 1]), ])
    colnames(y.o) <- c("time", "status")
    
    # uncensored patients
    I <- which(y.o$status==1)
    # patients at risk
    R <- lapply( which(y.o$status==1) , function(i) which( y.o$time >= y.o$time[i] ) )
    names(R) <- paste0("R", which(y.o$status==1))
    
    lambda.max <- rep(0, B)
    eps <- 0.01 # cf Yang et al (2013) the cocktail algorithm
#     for (i in 1:B){
#       lambda.max[[i]] <- 0.5*sum(apply(X[[i]], 2, function(c) norm(c, "2")))
#     }
#     a <- b <- NULL
#     for (i in 1:B){
#       b[[i]] <- apply(X.o[[i]], 2, function(v) sapply(names(R), function(r) {abs(sum(-v[R[[r]]] + 1/length(R[[r]]) * sum(v[R[[r]]]) ))} ))
# #       a[[i]] <- abs(b[[i]])
#       lambda.max[[i]] <- (1/n)*(1/alpha) * max(b[[i]])
# #       lambda.max[[i]] <- (1/n)*(1/alpha) * max(apply(X.o[[i]], 2, function(v) abs(sum(sapply(names(R), function(r) {-v[R[[r]]] + 1/length(R[[r]]) * sum(v[R[[r]]]) } )))))
#     }
### avec boucles for pour block 1
d <- NULL
for (b in 1:B){
  d[[b]] <- rep(0, ncol(X.o[[b]]))
  for (i in 1:length(names(R))){
    a <- matrix(0, nrow=length(R[[names(R)[[i]]]]), ncol=ncol(X.o[[b]]));
    a <- -X.o[[b]][R[[names(R)[i]]], ] + (1/length(R[[names(R)[i]]]) * sum(X.o[[b]][R[[names(R)[i]]], ]) )
    for (j in 1:ncol(a)){
      d[[b]][j] <- abs(sum(a[,j]))
    }
  }
  lambda.max[[b]] <- max(d[[b]])/(n*alpha)
}
    lambda.min <- eps*lambda.max
#     print(b)
#     print(a)
#     print(paste("lambda max : ", lambda.max))
#     print(paste("lambda min : ", lambda.min))
    # -1 because sequence starts at 0 (0:nb.lambda) and -1 to include huge lambda (<Inf)
    nb.lambda <- nl - 1 #- 1 # corresponds to m, in the paper
    l <- NULL
    # 0 to generate keep space for huge lambda
    j <- c(0, 0:nb.lambda) 
    for(i in 1:B){
      l[[i]]<- lambda.max[[i]]*(lambda.min[[i]]/lambda.max[[i]])^(j/nb.lambda)
      # introduce huge lambda at the beginning
#       l[[i]][[1]] <- 1e+4
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
