source("/home/philippe/github/multiblox/istacox_MapRedR/scripts/functions.R")

sparse.partial.loglik <- function(model, newdata, newy, lambda){
  # model : beta of the fitted model
  # newdata : data from de test set
  # newy : outcome from de the test set to compute I and R
  
  beta <- model$beta
  ### Sets of patient at risk at Ti
  x.o <- newdata[order(newy[, 1]), ]
  y.o <- as.data.frame(newy[order(newy[, 1]), ])
  
  # uncensored patients
  I <- which(y.o$status==1)
  # patients at risk
  R <- lapply( which(y.o$status==1) , function(i) which( y.o$time >= y.o$time[i] ) )
  names(R) <- paste0("R", which(y.o$status==1))
  
  pll <- mapply( function(i, j) sum(newdata[i, ]%*%beta - log(sum( exp(newdata[R[[sprintf("R%i", j)]], ]%*%beta) ))), I, R)
  alpha <- 0.5*lambda
  gamma <- 0.25*lambda
  spll <- pll - (alpha*norm.l1(beta)) - (gamma*norm.l2(beta))
  
  return(spll=spll)
}