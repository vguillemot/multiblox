#' Performance score
#' 
#' @param y.test
#' @param y.hat
#' @return performance score
#' @keywords internal
istacox.score <- function(y.test, y.hat){
    return(list(perf=y.hat))
}
