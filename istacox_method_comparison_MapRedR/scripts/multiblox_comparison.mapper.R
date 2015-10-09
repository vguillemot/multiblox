# Mapper for multiblox_comparison
# Calls multiblox.lambda.tune

args <- commandArgs(trailingOnly = TRUE)
print(args)

# common parameters from command line : fold_file, map_output_file, grid.file, k, i, inner_cv_method, inner_cv_metric, des, pathtoscript
# Rscript /home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/multiblox_comparison.mapper.R 
# 1) /home/philippe/Projets/these/multiblox/multiblox_methods_comparison_MapRedR_results/test/fold_data_1.Rdata 
# 2) /home/philippe/Projets/these/multiblox/multiblox_methods_comparison_MapRedR_results/test/hierarchical.map_res_1_1.Rdata 
# 3) /home/philippe/Projets/these/multiblox/multiblox_methods_comparison_MapRedR_results/test/grid.multiblox_comparison.Rdata
# 4) 1 
# 5) 1 
# 6) CV 
# 7) deviance 
# 8) hierarchical 
# 9) /home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/
fold_file <- args[1]
output_file <- args[2]
grid.file <- args[3]
inner_fold <-strtoi(args[4])
outer_fold <- strtoi(args[5])
inner_cv_method <- args[6]
inner_cv_metric <- args[7]
# custom parameters : model.custom.param.map
design <- args[8]
pathtoscript <- args[9]

adaptative <- F
fast <- F

print(paste("design : ", design, sep=""))
print(paste("adaptative : ", adaptative, sep=""))
print(paste("fast : ", fast, sep=""))
print(paste("pathtoscript : ", pathtoscript, sep=""))

# load lambda.grid
load(grid.file)
lambda.grid <- parameters.grid

if (is.na(pathtoscript)) { pathtoscript = "./" }
source(paste(pathtoscript, "multiblox.lambda.tune.R", sep=""))

load(fold_file)

if (design == "hierarchical") { 
  D <- matrix(0, ncol=B, nrow=B) 
} else if (design == "complete") { 
  D <- matrix(1, ncol=B, nrow=B) - diag(B)
}

#X, y, D, trainmat, i, outer_it, lambda.grid, scale=T, method="CV"
reswr <- multiblox.lambda.tune(X=X.train, y=y.train, D=D, trainmat, inner_fold, 
                              outer_fold, lambda.grid, scale=scale, method=inner_cv_method, metric=inner_cv_metric, adaptative=adaptative, fast=fast)

cur.model <- reswr[["res"]]
cur.pred.score <- as.matrix(reswr[["pred.score"]])
save(cur.pred.score, lambda.grid, cur.model, design, file=output_file)
cat("Writing", output_file, "\n")