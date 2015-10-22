##################################################
## A multiblock approach to Cox regression      ##
###           --  simulated data  --           ###
##################################################

rm(list=ls()) ; gc()

# Load libraries.
library(mvtnorm)
library(CMA)
library(MULTIBLOX)
library(rjson)

#########################################################################
###
###  reading JSON configuration file, parameters setting and data loading
###
#########################################################################
args <- commandArgs(trailingOnly = TRUE)

# Set working directory.
data_json_file <- args[1]
cv_json_file <- args[2]
data.cfg <- fromJSON(paste(readLines(data_json_file), collapse=""))
cv.cfg <- fromJSON(paste(readLines(cv_json_file), collapse=""))

pathtowritefile <- data.cfg$pathtowritefile
pathtoscript <- cv.cfg$pathtoscript
# source required scripts
# sources <- c("functions.R", 
#              "generate_learning_sets.R")
# for(s in sources) { source(paste(pathtoscript, s, sep="")) }

# CV parameters
nf <- cv.cfg$nf   # nb of folds for outer CV
nfi <- cv.cfg$nfi # nb of folds for inner CV
if (is.null(nfi)) { nfi <- 10 }
cv_method <- cv.cfg$cv_method
cv_metric <- cv.cfg$cv_metric
inner_cv_method <- cv.cfg$inner_cv_method
inner_cv_metric <- cv.cfg$inner_cv_metric

# Selects data source and reads corresponding parameters
source(data.cfg$data_source)
data <- read.data(cfg=data.cfg)
X <- data$X
y <- data$y
B <- length(X)
my_var_names <- data$my_var_names


set.seed(42)
inner_cv_seed <- 4257
scale <- T

#########################################################################
###
###  write chunks of data
###
#########################################################################
### set parameter for CV
N <- nrow(X[[1]])

### Generation of the training and testing sets with a method from package CMA.
### trainig sets are stratified according to the censoring variable y[, 1]
if (cv_method =="MCCV") { # with 0.8 train and 0.2 test
  trainmat.outer <- GenerateLearningsets(y=as.factor(y[, 2]),method="MCCV", niter=nf, ntrain=floor(N * 0.8), strat=TRUE)@learnmatrix
} else {
  trainmat.outer <- GenerateLearningsets(y=as.factor(y[, 2]),method="CV",fold=nf, strat=TRUE)@learnmatrix
}

# Number of inner folds equal number of subjects in the train set (outer CV).
if (inner_cv_method == "LOOCV") { nfi <- ncol(trainmat.outer) }

# store original data anf nfi in case of LOOCV
original_file = paste(pathtowritefile, "original_data.Rdata", sep="")
save(data, nfi, file=original_file)
print(original_file)

### block concatenation
X.concat <- as.data.frame(Reduce("cbind", X))
colnames(X.concat) <- unlist(my_var_names)

### Outer cross validation loop
for (i in 1:nrow(trainmat.outer)){
  ##1) get train and test for this iteration
  ind <- trainmat.outer[i,]
  
  # reindexing
  if (cv_method=="LOOCV") {
    X.test.concat <- t(as.matrix(X.concat[-ind,]))
  } else {
    X.test.concat <- X.concat[-ind, ]
  }
  X.test <- lapply(X, function(mm) { as.matrix(mm[-ind,])})
  X.train.concat <- X.concat[ind,]
  X.train <- lapply(X, function(mm) mm[ind, ])
  y.train <- y[ind,]
  y.test <- y[-ind,]
  # scaling
  if (scale) {
    # Comment: concat scaled in map_glmnet.R
    X.train <- lapply(X.train, function(mm) scale2(mm))
    scl_fun <- function(data, scaled) {
      scale(data, center = attr(scaled, "scaled:center"),
            scale = attr(scaled, "scaled:scale")) }
    X.test <- mapply(scl_fun, X.test, X.train, SIMPLIFY=FALSE)
  }
  # generate inner folds
  N.train = nrow(X.train[[1]])
  trainmat <- generate.learning.sets(N.train, y.train[, 2], method=inner_cv_method, nf.cv=nfi, inner_cv_seed=inner_cv_seed)
  
  # save data for the fold
  fold_file = paste(pathtowritefile, "fold_data_", i, ".Rdata", sep="")
  save(X.train, y.train, X.test, y.test, X.train.concat, X.test.concat, cv_method, cv_metric, 
       scale, trainmat, B, file=fold_file)
  print(fold_file)
}
