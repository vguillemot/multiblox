norm.l1 <- function(beta){
  return(sum(abs(beta)))
}

norm.l2.2 <- function(beta){
  return(sum(beta^2))
}
