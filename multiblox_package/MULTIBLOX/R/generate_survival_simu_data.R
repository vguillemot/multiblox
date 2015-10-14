#' Generate multi-block survival simulated data.
#' @param n total number of individuals.
#' @param B number of non-CGH blocks.
#' @param B.CGH number of CGH blocks.
#' @param p vector of number of variables per block.
#' @param n.expl vector of number of informative variables per block.
#' @param noise_corr vector of correlation noise parameters per block.
#' @param alpha ???
#' @param seed random seed.
#' @return Simulated data: X, y, names of the variables and hazard function.
#' @examples 
#' nf <- 10 # nb of folds for CV
#' n <- nf*nf*2 #observations
#' p <- c(40, 30) # c(40, 30, 50, 20, 30, 45) #nb of variables in each block
#' noise_corr <- c(0.8, 0.5) # correlation between variables within each block
#' B <- length(p) # number of blocks
#' rho <- 0.7 # correlation parameter !!! choose carefully the rho's so that sigmaComponents is positive definite !!!
#' alpha <- 0.1 or sqrt(2) # coefficient that moderates eta's part in the model
#' B.CGH <- 1 # nb of CGH simulated blocks
#' nb.expl <- 100 # nb of variables with information in each block
#' library(multiblog) # for idn function
generate_survival_simu_data <-
function(n, B=2, B.CGH=1, rho=rep(0.7, B-1), p=rep(100, B), n.expl=rep(100, B), 
         noise_corr=rep(0.5, B), alpha=0.1, seed=4257){ #}, pathtodata="./"){
    
  require(MULTIBLOX) # à ajouter au début de chaque fichier
    # library(mvtnorm) # for rmvnorm function
    # library(cghseg) # to simulate CGH data (non-normal distribution)
    # library(glmnet)
    
#     if (B.CGH!=0){
#       source(paste(pathtodata, "simulCGHprofiles.R", sep="")) # hacked from cghseg function to set the desired seed
#     }
#     source(paste(pathtodata, "rmvnorm_withGivenLatent.R", sep=""))# hacked from rmvnorm function to provide CGH latent variables
#     source(paste(pathtodata, "functions.R", sep=""))
    
    set.seed(seed)
    
    B.norm <- B - B.CGH
    p.CGH <- p[(B.norm+1):B]
    
    if (B.CGH!=0){
      ### CGH data simulations  
      simu.CGH <- NULL
      pc1 <- matrix(0, nrow=n, ncol=B.CGH)
      for (i in 1:B.CGH){
        simu.CGH[[i]] <- t(as.matrix(simulCGHprofiles(M=n, n=p.CGH[i], k.mean=8, SNR=5, lambda=1, seed=seed)$Y))# hacked to set the desired seed
        pc1[,i] <- prcomp(simu.CGH[[i]])$x[, 1]
      }
    } else {
      pc1 <- matrix(0, nrow=n, ncol=B.norm)
      for (i in 1:B.norm){
        pc1[,i] <- rnorm(n = n)
      }
    }
    
    ### norm data simulations
    if (length(noise_corr)!=length(p)) stop("Attention, chaque bloc doit
                                            avoir une corrélation pour le bruit !")
    
    my_block_names <- paste0("x", 1:B) # c("x1", "x2", "x3", "x4", "x5", "x6")
    my_var_names <- lapply(1:B, function(b) sprintf("blk%dvar%d", b, 1:p[b]))
    k <- 1
    if(length(rho)==1){
      sigmaComponents <- toeplitz(c(1,rep(rho,B-1)))# covariance matrix between components
    } else if(length(rho)!= B*(B-1)/2){
      stop("Attention, la longueur de rho est soit 1 soit (B*(B-1))/2 !")      
    } else {
      sigmaComponents <- matrix(0, nrow=B, ncol=B)
      #       if (B == 2) {
      #         sigmaComponents[1, 1] <- 1
      #       } else {
      for (i in 1:(B-1)){
        for (j in (i+1):B){
          #print(k)
          sigmaComponents[i, j] <- rho[k]
          k <- k+1
          if(j==B) {
            break
          }
        }
      }
      #       }
      #print(sigmaComponents)
      I <- idn(sigmaComponents)
      sigmaComponents <- sigmaComponents + t(sigmaComponents) + I
    }
    #print(sigmaComponents)
    etaks <- rmvnorm_withGivenLatent(n=n, latent=pc1, sigma=sigmaComponents) # latent component 1 per non-y block, k=block index
    sigmaz <- mapply(noise_corr, p, FUN=function(nc, q) toeplitz(rep(c(1,nc),c(1,q-1))), SIMPLIFY=FALSE)# covariance matrix of noise for each block
    beta <- NULL
    alphas <- NULL
    z <- x <- eta <- NULL
    for (b in 1:B){
      alphas[[b]] <- rep(c(alpha,0),c(n.expl[b], p[b]-n.expl[b])) # n.expl[b] variables contain information
      z[[b]] <- rmvnorm(n=n, mean=rep(0,p[[b]]), sigma=sigmaz[[b]]) # library(mvtnorm) # noise in block b
      x[[b]] <- etaks[,b]%*%t(alphas[[b]]) + z[[b]] # data in block b
    }
    
    if (B.CGH!=0){
      for (d in (B.norm+1):B){
        x[[d]] <- simu.CGH[[d-B.norm]]
      }
    }
    
    for (c in 1:B){
      beta[[c]] <- seq(from=-2*c, to=2*c, length.out=p[c]) + rnorm(p[c], 0, c)
      eta[[c]] <- x[[c]]%*%beta[[c]]
    }
    combi <- 0*rnorm(n, 0, 1)
    for (b in 1:B){
      print(eta[[b]])
      combi <- combi + eta[[b]]
    }
    # code inspired from glmnet
    ## ep = event proportion
    ep <- 0.5
    tcens=rbinom(n=n,prob=ep,size=1) # n values in binomial distribution = censoring indicator
#     y.o <- 1-tcens #on suppose que les observations sont ordonnées
#     
#     # uncensored patients
#     I <- which(y.o==1)
#     # patients at risk
#     R <- lapply( which(y.o==1) , function(i) return(c(i:n)) )
#     names(R) <- paste0("R", which(y.o==1))
    h0 <- 0.2
#     for(b in 1:B){
#       h0[[b]] <- mapply( function(i, j) sum(y.o[j] / sum( exp(x[[b]][R[[sprintf("R%i", i)]], ]%*%beta[[b]]))), I, R) # baseline hazard
#     }
#     
#     print(unlist(h0))
    hx=h0*exp(combi/n) # hazard function
#     print(hx)
    ty=rexp(n,hx) # n values in exponential distribution = survival time 
#     print(ty)
    y=data.frame(time=ty,status=1-tcens, stringsAsFactors = F) # y=Surv(ty,1-tcens) with library(survival)
    
    X=list()
    if(B==1){
      X[[1]] = x[[1]]
      names(X)=my_block_names[1:B]
      colnames(X[[1]]) <- my_var_names[[1]]
    }else{
      X=x
      names(X)=my_block_names[1:B]
      for (i in 1:B){
        colnames(X[[i]]) <- my_var_names[[i]]
      }
    }
    
    
    return (list(X=X, y=as.matrix(y), my_var_names=my_var_names, h0=h0))
  }
