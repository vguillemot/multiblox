#' Performs block relaxation (ex-ALS Alternating Least Squares)
#' 
#' @param x is a list of B matrices (blocks)
#' @param I list of non censored individuals, ranked by event time.
#' @param R list of sets of individuals at risk at each non censored ranked time
#' @param D is the matrix of design storing the links among blocks
#' @param lambda a tuple of L1-norm shrinkage parameter
#' @param eps is the required precision.
#' @param max.iter is the maximal number of block relaxation iterations.
#' @param beta_init is the initial value of the weight vector (used mainly for warm restarts).
#' @param fast is a boolean used to specify if FISTA (TRUE) is used instead of ISTA (FALSE, default).
#' @param ada is a boolean used to specify is the step must be chosen at each iteration (TRUE) or not (FALSE, default).
#' @return beta a list of multiblox coefficients, convergence speed of convergence, niter nb of iterations
#' @keywords internal
relax_multiblox <-
function(x, I, R, D, lambda = 0, eps = 0.001, max.iter = 10000, beta.init = NULL, fast=fast, ada=ada){
  ### Block relaxation for multiblox
  
    # x is a list of B (p_k by n) matrices of predictors
    # I vector of indexes of uncensored patients, 
    # R list of vector of indexes of patients at risk for each ordered time Ti, 
    # D is a B by B binary design matrix describing connections between blocks
    # lambda.opt is a vector of shrinkage parameters, one for each block
    # eps is a real corresponding to the tolerance for convergence
    # max.iter is an integer corresponding to the maximum nb of ALS iterations
    # beta.init is a p by 1 vector of the initial values of beta
    
#     source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/istacox.R")
#     source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/link.R")
    
    B <- length(x) #number of blocks
    n <- nrow(x[[1]])
    p <- sapply(x, ncol) #number of covariates in each block
    
    beta <- beta_new <- e <- S <- eta <- cvg <- link <- lk <- NULL
    
    ### Initialization step
    iter <- 1
    cvg[[iter]] <- 10
    beta <- list()
    
    if (!is.null(beta.init)){
      beta <- beta_new <- beta.init
    }else{
      f <- function(x) x/1e2
      for (c in 1:B){
        beta[[c]] <- rep(0, p[c])
      }
      beta_new <- beta
    }
    iter.inner <- 0
    div.inner <- 0
    max.iter.inner <- 1000
    consec_max.iter <- 0
    beta_new <- beta
    crit <- pen <- mb_crit <- 0
    
    print(D)
    
    for (iter in 1:max.iter){
      print(paste("relax iter : ", iter, sep=""))
      for (b in 1:B) {
#         print(paste("block : ", b))
#         print(dim(x[[b]]))
#         print(dim(beta[[b]]))
        link <- link(x, D, b, beta)
        istacox_res <- istacox(X=x[[b]], I=I, R=R, alpha=0.5*lambda[[b]], kmax=1000, epsilon=1e-4, 
                           fast=fast, ada=ada, link=link, beta_init=beta[[b]])
        iter.inner <- iter.inner + istacox_res$k
        print(paste("istacox iter : ", iter.inner, sep=""))
        if (max.iter.inner == istacox_res$k) {
          div.inner <- div.inner + 1
          consec_max.iter <- consec_max.iter + 1
        } else {
          consec_max.iter <- 0
        }
#               beta_new <- lapply(NR_step$beta, beta_norm2)
        beta_new[[b]] <- istacox_res$beta
      }
print("un tour de block relaxation")
      d <- mapply("-", beta, beta_new, SIMPLIFY=FALSE)
      e[[iter]] <- sapply(d, base::norm, "f")
      if(max(unlist(e))<eps) break
      beta <- beta_new
### Calcul du critère
      for (c in 1:B){
        pll <- sum(mapply( function(i) {x[[c]][i, ]%*%beta[[c]] - log(sum( exp(as.matrix(x[[c]][R[[sprintf("R%i", i)]], ])%*%beta[[c]]) ))}, I))
        pen <- pen + lambda[[c]]*(norm.l1(beta[[c]])+0.5*norm.l2.2(beta[[c]]))
        for (d in setdiff(1:B, c)){
          lk[[c]] <- sum(D[c, d]*t(x[[c]]%*%beta[[c]])%*%x[[d]]%*%beta[[d]])
        }
        crit <- crit+pll
      }
      mb_crit <- - crit - sum(lk) + pen
# print(dim(crit))
# print(dim(link))
# print(dim(pen))
# print(dim(mb_crit))
print(paste("critere : ", mb_crit, sep=""))
      if (consec_max.iter >= (2 * B)) break # all the blocks diverge 2 times
    }
    #eta <- mapply(x, beta, FUN=function(a, b) a%*%b, SIMPLIFY=FALSE) ### bug 3 blocs
    print(paste("Block relaxation iter (inner - div): ", iter , "(", iter.inner, "-", div.inner, ")" ))
    return(list(beta = beta_new, convergence = unlist(e), niter=iter, crit=crit))
}
