############################################################################################
##   A multiblock approach to coxnet regression for both designs hierarchical and complete  ###
###     Comparison with coxnet for ridge, lasso and elastic net shrinkage methods         ##
###                                   --  simulated data  --                             ###
############################################################################################
### call
# Rscript --vanilla main_generate_survival_data_chunks.R data.local.survival_2B_CGH_n64.json cv.scheme.local.survival_CGH_2B.json
# Rscript --vanilla main_multiblox_comparison.R data.local.survival_2B_CGH_n64.json method.multiblox.comparison.json | grep Rscript > cmd.comparison.n64.log
# python create_multiblox_workflow.py cmd.comparison.n64.log multiblox_comparison.n64.wf
# soma_workflow_gui&
# open multiblox_comparison.n64.wf
# submit
###################################################

rm(list=ls()) ; gc()

# Load libraries.
library(mvtnorm)
library(CMA)
library(MULTIBLOX)
library(rjson)
library(glmnet)

#########################################################################
###
###  reading JSON configuration file, parameters setting and data loading
###
#########################################################################
args <- commandArgs(trailingOnly = TRUE)

# Set working directory.
data_json_file <- args[1]
method_json_file <- args[2] # with any design because it will perform both
data.cfg <- fromJSON(paste(readLines(data_json_file), collapse=""))
method.cfg <- fromJSON(paste(readLines(method_json_file), collapse=""))

pathtowritefile <- data.cfg$pathtowritefile
pathtoscript <- method.cfg$pathtoscript
# source required scripts
# sources <- c("multiblox.lambda.tune.R", "functions.R", "make.lambda.grid.R", "generate_learning_sets.R")
# for(s in sources) { source(paste(pathtoscript, s, sep="")) }

# method
model.name <- method.cfg$model.name
model.mapper <- method.cfg$model.mapper
model.inner.reducer <- method.cfg$model.inner.reducer
print("voici les parametres : ")
print(model.name)
print(model.mapper)
print(model.inner.reducer)

# CV parameters
nf <- method.cfg$nf   # nb of folds for outer CV
nfi <- method.cfg$nfi # nb of folds for inner CV
if (is.null(nfi)) { nfi <- 10 }
cv_method <- method.cfg$cv_method
cv_metric <- method.cfg$cv_metric
inner_cv_method <- method.cfg$inner_cv_method
inner_cv_metric <- method.cfg$inner_cv_metric

# method parameters
method.map.param <- method.cfg$model.custom.param.map # fast adaptative pathtoscript
method.ired.param <- method.cfg$model.custom.param.ired # fast adaptative pathtoscript
print(method.map.param)
print(method.ired.param)

coxnet_res_file <- paste(pathtowritefile, "coxnet_red_res_",sep="")
scale <- T

# data
# Selects data source and reads corresponding parameters
source(data.cfg$data_source)
data <- read.data(cfg=data.cfg)
X <- data$X
y <- data$y
B <- length(X)
my_var_names <- data$my_var_names

# construct the parameters grid
parameters.grid <- make.lambda.grid(X, path="smart")
# print(lambda.grid)
grid.file <- paste(pathtowritefile, "grid.",  model.name, ".Rdata", sep="")
save(parameters.grid, file=grid.file)

set.seed(42)
designs <- c("hierarchical", "complete")
do_coxnet <- TRUE
inner_cv_seed <- 4257

#########################################################################
###
###  prediction performances assessment
###
#########################################################################
### set parameter for CV
N <- nrow(X[[1]])

### Generation of the training and testing sets with a method from package CMA.
### TODO: generate.learning.sets
if (cv_method =="MCCV") {
  trainmat.outer <- GenerateLearningsets(y=as.factor(y[, 2]),method="MCCV", niter=nf, ntrain=floor(N * 0.8), strat=TRUE)@learnmatrix
} else {
  trainmat.outer <- GenerateLearningsets(y=as.factor(y[, 2]),method="CV",fold=nf, strat=TRUE)@learnmatrix
}

# Number of inner folds equal number of subjects in the train set (outer CV).
if (inner_cv_method == "LOOCV") { nfi <- ncol(trainmat.outer) }

### block concatenation
X.concat <- as.data.frame(Reduce("cbind", X))
colnames(X.concat) <- unlist(my_var_names)

# glmnet call for coxnet !!
if (do_coxnet) {
  for (i in 1:nrow(trainmat.outer)){
    cat(paste("Rscript",  paste(pathtoscript, "coxnet.mapper.R", sep=""), paste(pathtowritefile, "fold_data_",  i, ".Rdata", sep=""), 
            paste(coxnet_res_file, i, ".Rdata", sep=""), nf, cv_method, cv_metric, "\n", sep=" ")) 
  }
#   cat(paste("Rscript",  paste(pathtoscript, "coxnet.outer.reducer.R", sep=""), 
#              coxnet_res_file, nf, paste(pathtowritefile, "coxnet_res.Rdata", sep=""), "\n",sep=" ")) 
}

### Cross validation loop
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
    # Comment: concat scaled in coxnet.mapper.R
    X.train <- lapply(X.train, function(mm) scale2(mm))
    scl_fun <- function(data, scaled) {
      scale(data, center = attr(scaled, "scaled:center"),
            scale = attr(scaled, "scaled:scale")) }
    X.test <- mapply(scl_fun, X.test, X.train, SIMPLIFY=FALSE)
  }
  # generate inner folds STRATIFICATON ACCORDING TO STATUS
  N.train = nrow(X.train[[1]])
  trainmat <- generate.learning.sets(N.train, y.train[, 2], method=inner_cv_method, nf.cv=nfi, inner_cv_seed=inner_cv_seed)
  
  # save data for the fold
  colnames(y.train) <- colnames(y.test) <- c("time", "status")
  fold_file = paste(pathtowritefile, "fold_data_", i, ".Rdata", sep="")
  save(X.train, y.train, X.test, y.test, X.train.concat, X.test.concat, cv_method, cv_metric, 
       scale, trainmat, parameters.grid, B, file=fold_file)
  
  
  outer_red_pattern <- paste(pathtowritefile, model.name, ".red_res_" , sep="")
  
  # generate multiblox commands for inner folds 
  # no inner fold for coxnet because it has its own model selection procedure
  for (des in designs) {    ### pour chaque design hierarchical and complete
    print(des)
    ### Cross validation loop
      # data for the fold
      outer_res_file = paste(outer_red_pattern, des, i, ".Rdata", sep="")
        # generate inner CV commands
        map_file_pattern = paste(pathtowritefile, des,".map_res_", i, "_", sep="")
        for (k in 1:nfi){
          map_output_file = paste(map_file_pattern, k, ".Rdata", sep="")
          ### commandes de mapper multiblox
          cat(paste("Rscript", paste(pathtoscript, model.mapper, sep=""), fold_file, 
                    map_output_file, grid.file, k, i,
                    inner_cv_method, inner_cv_metric, des, method.map.param, "\n",sep=" "))
        }
        ### commandes inner reducer multiblox
        cat(paste("Rscript", paste(pathtoscript, model.inner.reducer, sep=""), fold_file, map_file_pattern, outer_res_file,
                  i, nfi, des, method.ired.param, "\n",sep=" "))
  }
}

## final reducer call
red_pattern <- paste(pathtowritefile, "multiblox.red_res_" , sep="")
res_file <- paste(pathtowritefile, "final_res.Rdata" , sep="") 
### commandes outer reducer qui calcule les stats pour chaque méthode.
cat(paste("Rscript", paste(pathtoscript, "multiblox_comparison.outer.reducer.R", sep=""), coxnet_res_file, red_pattern, nf, 
          res_file, "\n",sep=" "))

