#' Centers and scales a matrix A
#' 
#' @param A is a matrix of observations
#' @param center is a boolean indicating if centering should be performed
#' @param scale is a boolean indicating if scaling should be performed
#' @return the new matrix A, centered and/or scaled, with scaling factor as an attribute
#' @keywords internal
scale2 <-
function (A, center = TRUE, scale = TRUE) 
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
