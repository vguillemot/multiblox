#' Computes a grid of parameter
#' 
#' @param data is the matrix of observations
#' @return a grid of lambdas, sparsity parameter
#' @keywords internal
make.grid <-
function(data){
  return (make.lambda.grid(data$X, path="naive"))   
}
