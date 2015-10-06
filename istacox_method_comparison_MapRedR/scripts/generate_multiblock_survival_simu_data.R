generate_multiblock_survival_simu_data <- function(n=n, p=p, B=B, p_frac=rep(0.25, B), corr=0.5){
  require(survJamda)
  
  set.seed(4256)
  ds <- X <- my_var_names <- NULL
  pp <- sum(p)
  ff <- sum(p*p_frac)
    
  ds <- generate.survival.data(gene.nb = ff, tot.genes = pp, correlation = corr,
                                 sample.nb = n, beta.init = 1.5, shape = 1, scale = 1)
  
  ### divide the super block into B blocks and allocate the variables among the blocks
  ### meaningful variables are at the beginning of the each block.
  var <- rep(1:B, p*p_frac)
  no_var <- rep(1:B, p-p*p_frac)
  
  for (b in 1:B){
    my_var <- ds$ds1[, which(var==b)]
    my_no_var <- ds$ds1[, which(no_var==b)]
    X[[b]] <- cbind(my_var, my_no_var)
    print(dim(X[[b]]))
    print(p[b])
    my_var_names[[b]] <- sprintf("blk%dvar%d", b, 1:p[b])
    colnames(X[[b]]) <- my_var_names[[b]]
  }
  
  return(list(X=X, y=as.matrix(cbind(ds$T, ds$censor)),  my_var_names=my_var_names, h0=1))
}


read.data <-
  function(cfg){
    n <- cfg$simu$n
    p <- cfg$simu$p
    p.frac <- cfg$simu$p.frac
    B <- cfg$simu$B
    if(length(cfg$simu) == 0) {
      print("bad data configuration")
      quit("no")
    }
    # basic configuration checking
    if ((is.null(n)) + (is.null(p)) + (is.null(p.frac))) {
      print("Missing or bad parameters to generate simulated data")
      quit("no")
    }
    res <- generate_multiblock_survival_simu_data(n=n, p=p, B=B, p_frac=p.frac, corr=0.5)
    return (res)
  }
