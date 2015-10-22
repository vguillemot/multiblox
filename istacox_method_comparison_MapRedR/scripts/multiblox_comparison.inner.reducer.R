# reducer for model selection in multiblox
# source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/istacox.R")
# source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/istacox.predict.R")
# source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/istacox.score.R")
# source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/istacox.R")
# source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/istacox.predict.R")
# source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/istacox.score.R")

library("MULTIBLOX")
args <- commandArgs(trailingOnly = TRUE)

# common parameters from command line : Rscript model.inner.reducer fold_file map_file_pattern outer_res_file i nfi cv_metric design
# 1) /home/philippe/Projets/these/multiblox/multiblox_methods_comparison_MapRedR_results/test/fold_data_1.Rdata 
# 2) /home/philippe/Projets/these/multiblox/multiblox_methods_comparison_MapRedR_results/test/hierarchical.map_res_1_ 
# 3) /home/philippe/Projets/these/multiblox/multiblox_methods_comparison_MapRedR_results/test/multiblox_comparison.red_res_hierarchical1.Rdata
# 4) 1 
# 5) 5 
# 6) hierarchical
input_file_fold <- args[1]
input_file_pattern <- args[2]
outputfile <- args[3]
outer_fold <- strtoi(args[4])
n_inner_folds <- strtoi(args[5])
design <- args[6]
# custom parameters
fast <- args[7]
adaptative <- args[8]
pathtoscript <- args[9]

print(adaptative)
print(fast)

load(input_file_fold) # data, cv_method, cv_metric, ...

### Beta initialization
#beta0 <- matrix(0, nrow=(ncol(X.train)), ncol=1)

### initialization
load(paste(input_file_pattern, 1,".Rdata", sep="")) # cur.pred.score (spll) 1 x nb(lambda) cur.model(list[[nb(lambda)]]) lambda.grid
n_lambdas = length(lambda.grid)
pred.score <- matrix(NA, nrow=n_lambdas, ncol=n_inner_folds) 
models <- NULL

### il faut calculer la vraisemblance partielle cross-validée
### pour cela, il nous faut le training set de la boucle externe de CV
### CALCUL DE LA VRAISEMBLANCE CROSS-VALIDEE

# lambda.opt <- unlist(lambda.grid[which.max(rowSums(pred.score, na.rm=T)), ])
# cat("Lambda optimal: ", lambda.opt, "\n")
#ret <- list(opt = lambda.opt, pred.score=pred.score, models=models, lambda.grid=lambda.grid)

### Sets of patient at risk at Ti
if(B==1){
  x.o <- X.train[[1]][order(y.train[, 1]), ]
}else{
  x.o <- lapply(X.train, function(l) l[order(y.train[, 1]), ])
}
y.o <- as.data.frame(y.train[order(y.train[, 1]), ])

print(dim(x.o))

# uncensored patients
I.train <- which(y.o$status==1)
# patients at risk
R.train <- lapply( which(y.o$status==1) , function(i) which( y.o$time >= y.o$time[i] ) )
names(R.train) <- paste0("R", which(y.o$status==1))

### gather all results
### dans le cas de multiblox, istacox.lambda.tune ne calcule que la vraisemblance partielle du training set (pred.score)
### et renvoie en plus le modèle
CV <- NULL
total.CV <- matrix(NA, nrow=nrow(lambda.grid), ncol=n_inner_folds)
for(i in 1:n_inner_folds) {
  inputfile = paste(input_file_pattern, i,".Rdata", sep="")
  cat("DO:", inputfile, "\n")
  load(inputfile)
  for(l in 1:nrow(lambda.grid)){
    #   models[[i]] <- cur.model # pour tous les lambdas de la grille
    #   pred.score[, i] <- cur.pred.score # pour tous les lambdas de la grille
    # print(lapply((lapply(x.o, function(l) l[R.train[[sprintf("R%i", i)]], ])), dim))
#     CV[l,i] <- sum(lapply(x.o, function(m) 
#       {mapply(function(i) 
#         {matrix(m[i, ], nrow=1) %*% matrix(cur.model[[l]][["beta"]], ncol=1) - log( sum( exp(m[R.train[[sprintf("R%i", i)]], ] %*% cur.model[[l]][["beta"]])))}, I.train)}))
#     - cur.pred.score[l,1]
    # print(lambda.grid[l, ])
    for (b in 1:B){
      CV[[b]] <- sum(mapply(function(i) {x.o[[b]][i,] %*% cur.model[[l]][["beta"]][[b]] - log( sum( exp(x.o[[b]][R.train[[sprintf("R%i", i)]], ] %*% cur.model[[l]][["beta"]][[b]])))}, I.train)) - cur.pred.score[l,1]
    }
    total.CV[l, i] <- sum(unlist(CV))
    #pll <- sum(mapply( function(i) newdata[i, ]%*%beta - log(sum( exp(newdata[R[[sprintf("R%i", i)]], ]%*%beta) )), I))
  }
}
# print(total.CV)

if (design == "hierarchical") { 
  D <- matrix(0, ncol=B, nrow=B) 
} else if (design == "complete") { 
  D <- matrix(1, ncol=B, nrow=B) - diag(B)
}

### model selection one lambda parameter for each block
lambda.opt <- lambda.grid[which.max(rowSums(total.CV, na.rm=T)),]
print(lambda.opt)
### refit
model.refit <- relax_multiblox(x=x.o, I=I.train, R=R.train, D, lambda=lambda.opt, max.iter=1000, eps=1e-4, beta.init=NULL,
                       fast=as.logical(fast), ada=as.logical(adaptative))
cur.beta.train <- model.refit[["beta"]]

### Pour l'evaluation du modele, il faudra calculer :
### 1) la déviance
### 2) le pronostic index
print(cv_metric)
##3) PREDICT 
dev <- istacox.predict(model=cur.beta.train, x=X.test, y=y.test, D=D, lambda=lambda.opt, type="deviance")
pi <- istacox.predict(model=cur.beta.train, x=X.test, y=y.test, D=D, lambda=lambda.opt, type="pi")
model_dev <- dev[["est"]]
model_pi <- pi[["est"]]

##4) SCORE: get and accumulate the score
dev.score <- istacox.score(y.test=model_dev, y.hat=model_dev)[["perf"]]
pi.score <- istacox.score(y.test=model_pi, y.hat=model_pi)[["perf"]]

save(dev.score, pi.score, lambda.opt, model.refit, design, file=outputfile)
cat("Writing", outputfile, "\n")