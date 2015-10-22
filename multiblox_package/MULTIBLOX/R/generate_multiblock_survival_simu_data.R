#' Generates multiblock simulated survival data
#' 
#' @param n number of observations.
#' @param p total number of variables.
#' @param B number of blocks.
#' @param p_frac fraction of variables with non-zero coefficients in the Cox model. 
#' @param corr correlation among variables.
#' @return X matrix of simulated data, y data-frame of survival data (time and status), my_var_names names of variables, h0 baseline hazard ratio

generate_multiblock_survival_simu_data <-
function(n=n, p=p, B=B, p_frac=rep(0.25, B), corr=0.5){
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
