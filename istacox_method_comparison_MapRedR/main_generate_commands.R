##################################################
## A multiblock approach to logistic regression ##
###           --  simulated data  --           ###
##################################################

rm(list=ls()) ; gc()

# Load libraries.
library(mvtnorm)
library(CMA)
#library(multiblog)
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
method_json_file <- args[3]
data.cfg <- fromJSON(paste(readLines(data_json_file), collapse=""))
cv.cfg <- fromJSON(paste(readLines(cv_json_file), collapse=""))
method.cfg <- fromJSON(paste(readLines(method_json_file), collapse=""))
# @TODO: better method file checking!
if(length(method.cfg) == 0) {
  print("bad method configuration")
  quit("no")
}

model.custom.param.map <- method.cfg$model.custom.param.map

pathtowritefile <- data.cfg$pathtowritefile
pathtoscript <- cv.cfg$pathtoscript
# source required scripts
sources <- c("istacox.lambda.tune.R", 
             "functions.R")
for(s in sources) { source(paste(pathtoscript, s, sep="")) }

# reload original data
load(paste(pathtowritefile, "original_data.Rdata", sep=""))

model.selection <- method.cfg$model.selection
model.name <- method.cfg$model.name

# construct parameters grid for this method
parameters.grid <- NULL
if(model.selection) {
  model.custom.param.ired <- method.cfg$model.custom.param.ired
  make.grid.func <- method.cfg$make.grid.func
  if(is.null(make.grid.func)) {
    print("bad method configuration")
    quit("no")
  } else {
    source(paste(pathtoscript, make.grid.func, sep=""))
    parameters.grid <- make.grid(data)
    # save grid
    grid.file <- paste(pathtowritefile, "grid.",  model.name, ".Rdata", sep="")
    if (file.exists(grid.file)) {
      warning(paste("Parameters grid file already exists. Change model.name or delete existing file.",
                    model.name, grid.file, sep=" "), immediate.=T)
    }
    save(parameters.grid, file=grid.file)
  }
}

# R scripts for MapReduce
model.mapper <- paste(pathtoscript, method.cfg$model.mapper, sep="")
model.outer.reducer <- paste(pathtoscript, method.cfg$model.outer.reducer, sep="")
if (model.selection) {
  model.inner.reducer <- paste(pathtoscript, method.cfg$model.inner.reducer, sep="")
}

# CV parameters
nf <- cv.cfg$nf   # nb of folds for outer CV

cv_method <- cv.cfg$cv_method
cv_metric <- cv.cfg$cv_metric
inner_cv_method <- cv.cfg$inner_cv_method
inner_cv_metric <- cv.cfg$inner_cv_metric
if (inner_cv_method == "LOOCV") { 
  # inherited from original_data.Rdata
  print(paste("nfi:", nfi))
  if (is.null(nfi)) { 
    print("Unable to retrieve nfi (number of inner folds)")
    quit("no")
  }  
} else {
  nfi <- cv.cfg$nfi # nb of folds for inner CV
}

#########################################################################
###
###  Generate command lines
###
#########################################################################
outer_red_pattern <- paste(pathtowritefile, model.name, ".red_res_" , sep="")
final_res_file <- paste(pathtowritefile, model.name, ".final_res.Rdata" , sep="") 
### Cross validation loop
for (i in 1:nf){
  # data for the fold
  fold_file = paste(pathtowritefile, "fold_data_", i, ".Rdata", sep="")
  outer_res_file = paste(outer_red_pattern, i, ".Rdata", sep="")
  if (model.selection) {
    # generate inner CV commands
    map_file_pattern = paste(pathtowritefile, model.name,".map_res_", i, "_", sep="")
    for (k in 1:nfi){
      map_output_file = paste(map_file_pattern, k, ".Rdata", sep="")
      cat(paste("Rscript", model.mapper, fold_file, map_output_file, grid.file, k, i,
                inner_cv_method, inner_cv_metric, model.custom.param.map, "\n",sep=" "))
    }
    cat(paste("Rscript",  model.inner.reducer, fold_file, map_file_pattern, outer_res_file,
              i, nfi, cv_metric, model.custom.param.map, "\n",sep=" "))
  } else {
    # write mapper command (without model selection)
    cat(paste("Rscript", model.mapper, fold_file, outer_res_file, i, model.custom.param.map, 
              "\n",sep=" "))
  }
}
## final reducer command
cat(paste("Rscript", model.outer.reducer, outer_red_pattern, nf, final_res_file, 
          "\n",sep=" "))

