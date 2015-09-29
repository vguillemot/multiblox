# Mapper for istacox
# Calls istacox.lambda.tune

args <- commandArgs(trailingOnly = TRUE)
print(args)

# common parameters from cross_validation_scheme.json and istacox_conf.json
fold_file <- args[1]
output_file <- args[2]
grid.file <- args[3]
inner_fold <-strtoi(args[4])
outer_fold <- strtoi(args[5])
method <- args[6]
cv_metric <- args[7]

# custom parameters
adaptative <- args[8]
fast <- args[9]
pathtoscript <- args[10]

print(paste("adaptative : ", adaptative, sep=""))
print(paste("fast : ", fast, sep=""))
print(paste("pathtoscript : ", pathtoscript, sep=""))

# load lambda.grid
load(grid.file)
lambda.grid <- parameters.grid

if (is.na(pathtoscript)) { pathtoscript = "./" }
source(paste(pathtoscript, "istacox.lambda.tune.R", sep=""))

load(fold_file)

#X, y, trainmat, i, outer_it, lambda.grid, scale=T, method="CV"
reswr <- istacox.lambda.tune(X=X.train, y=y.train, trainmat, inner_fold, 
                              outer_fold, lambda.grid, scale=scale, method=method, metric=cv_metric, adaptative=adaptative, fast=fast)
cur.model <- reswr[["res"]]
cur.pred.score <- as.matrix(reswr[["pred.score"]])
save(cur.pred.score, lambda.grid, cur.model, file=output_file)
cat("Writing", output_file, "\n")