

generate.learning.sets <-
function(N, y, method="MCCV", nf.cv=10, inner_cv_seed=4257) {
    
    library(CMA)
    set.seed(inner_cv_seed)
    if (nf.cv == N) {
      trainmat <- GenerateLearningsets(N, y, method = "LOOCV")@learnmatrix
      attr(trainmat, "method") <- "LOOCV"
    } else {
      if (method == "CV") {
        trainmat <- GenerateLearningsets(N, y, method = "CV", fold = nf.cv, strat=TRUE)@learnmatrix
        attr(trainmat, "method") <- "CV"
      } else {
        trainmat <- GenerateLearningsets(N, y, method="MCCV", niter=nf.cv, ntrain=floor(N * 0.8) , strat=TRUE)@learnmatrix
        attr(trainmat, "method") <- "MCCV"
      }
    }
  return(trainmat)
}
