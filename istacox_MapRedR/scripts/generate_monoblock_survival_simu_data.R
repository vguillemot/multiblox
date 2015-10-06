generate_monoblock_survival_simu_data <- function(n=n, p=p, p_frac=p.frac, corr=0.3){
  require(survJamda)
  
  set.seed(4256)
  ds <- generate.survival.data(gene.nb = p*p_frac, tot.genes = p, correlation = corr,
                         sample.nb = n, beta.init = 0.5, shape = 1, scale = 1)
  
  X=list()
  X[[1]] = ds$ds1
  my_var_names <- sprintf("blk1var%d", 1:p)
  colnames(X[[1]]) <- my_var_names
  
  return(list(X=X, y=as.matrix(cbind(ds$T, ds$censor)),  my_var_names=my_var_names, h0=1))
}


read.data <-
  function(cfg){
    n <- cfg$simu$n
    p <- cfg$simu$p
    p.frac <- cfg$simu$p.frac
    if(length(cfg$simu) == 0) {
      print("bad data configuration")
      quit("no")
    }
    # basic configuration checking
    if ((is.null(n)) + (is.null(p)) + (is.null(p.frac))) {
      print("Missing or bad parameters to generate simulated data")
      quit("no")
    }
    res <- generate_monoblock_survival_simu_data(n=n, p=p, p_frac=p.frac)
    return (res)
  }
