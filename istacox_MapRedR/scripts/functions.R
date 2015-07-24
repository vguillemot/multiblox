norm.l1 <- function(beta){
  return(sum(abs(beta)))
}

norm.l2.2 <- function(beta){
  return(sum(beta^2))
}

#' scale2() is function whose default method centers and scales the columns of a numeric matrix.
#' @param A  a numeric matrix
#' @param center  a logical value. If center = TRUE, each column is transformed to have zero mean (default: TRUE)
#' @param scale a logical value. If scale = TRUE, each column is transformed to have unit variance Value (default = TRUE)
#' @return list
#' {
#' \item{A}{The centered and/or scaled matrix. The centering and scaling values (if any) are returned as attributes "scaled:center" and "scaled:scale"}
#' } 
#'@export scale2


#INPUT
# A:      a numeric matrix
# center: a logical value. If center = TRUE, each column is transformed to have zero mean
# scale:  a logical value. If scale = TRUE, each column is transformed to have unit variance Value

#OUTPUT:
#A: The centered and/or scaled matrix. 
#   The centering and scaling values (if any) are returned as attributes "scaled:center" and "scaled:scale"


scale2<-function (A, center = TRUE, scale = TRUE) 
{
  if (center == TRUE & scale == TRUE) {
    A = scale(A, center = TRUE, scale = FALSE)
    std = apply(A, 2, sd)
    if (any(std==0)) {
      sprintf("there were %d constant variables",sum(std==0))
      std[std==0]=1
    }
    A = A/matrix(rep(std, nrow(A)), nrow(A), ncol(A), byrow = TRUE)
    attr(A, "scaled:scale") = std
    return(A)
  }
  if (center == TRUE & scale == FALSE) {
    A = scale(A, center = TRUE, scale = FALSE)
    return(A)
  }
  if (center == FALSE & scale == TRUE) {
    std = apply(A, 2, sd)
    A = A/matrix(rep(std, nrow(A)), nrow(A), ncol(A), byrow = TRUE)
    attr(A, "scaled:scale") = std
    return(A)
  }
}
