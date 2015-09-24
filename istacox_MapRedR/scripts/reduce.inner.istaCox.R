# reducer for model selection in multiblox
source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/istacox.R")
source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/istacox.predict.R")
source("/home/cathy/git_repo/multiblox/istacox_MapRedR/scripts/istacox.score.R")

args <- commandArgs(trailingOnly = TRUE)

# common parameters
input_file_fold <- args[1]
input_file_pattern <- args[2]
outputfile <- args[3]
outer_fold <- strtoi(args[4])
n_inner_folds <- strtoi(args[5])
# custom parameters
adaptative <- args[6]
fast <- args[7]

load(input_file_fold) # data, cv_method, cv_metric, ...

### Beta initialization
#beta0 <- matrix(0, nrow=(ncol(X.train)), ncol=1)

### initialization
load(paste(input_file_pattern, 1,".Rdata", sep=""))
n_lambdas = length(lambda.grid)
pred.score <- matrix(NA, nrow=n_lambdas, ncol=n_inner_folds) 
models <- NULL

### il faut calculer la vraisemblance partielle cross-validée
### pour cela, il nous faut le training set de la boucle externe de CV
### CALCUL DE LA VRAISEMBLANCE CROSS-VALIDEE

# lambda.opt <- unlist(lambda.grid[which.max(rowSums(pred.score, na.rm=T)), ])
# cat("Lambda optimal: ", lambda.opt, "\n")
ret <- list(opt = lambda.opt, pred.score=pred.score, models=models, lambda.grid=lambda.grid)

### Sets of patient at risk at Ti
if(B==1){
  x.o <- X.train[[1]][order(y.train[, 1]), ]
}else{
  x.o <- X.train[order(y.train[, 1]), ]
}
y.o <- as.data.frame(y.train[order(y.train[, 1]), ])

# uncensored patients
I.train <- which(y.o$status==1)
# patients at risk
R.train <- lapply( which(y.o$status==1) , function(i) which( y.o$time >= y.o$time[i] ) )
names(R.train) <- paste0("R", which(y.o$status==1))

### gather all results
### dans le cas de multiblox, istacox.lambda.tune ne calcule que la vraisemblance partielle du training set (pred.score)
### et renvoie en plus le modèle
CV <- NULL
for(i in 1:n_inner_folds) {
  inputfile = paste(input_file_pattern, i,".Rdata", sep="")
  cat("DO:", inputfile, "\n")
  load(inputfile)
  for(l in 1:length(lambda.grid)){
    #   models[[i]] <- cur.model # pour tous les lambdas de la grille
    #   pred.score[, i] <- cur.pred.score # pour tous les lambdas de la grille
    CV[l,i] <- sum(mapply(function(i){t(x.o[i, ]) %*% cur.model[[l]]$beta - log( sum( exp(x.o[R[[sprintf("R%i", i)]], ] %*% cur.model[[l]]$beta)))}, I)) - cur.pred.score[[l]]
  }

}
print(CV)

### choix du modele
lambda.opt <- unlist(lambda.grid[which.max(rowSums(CV, na.rm=T))])

### refit
model.refit <- istacox(X=x.o, I=I.train, R=R.train, alpha=0.5*lambda.opt, gamma=0.25*lambda.opt, kmax=1000, epsilon=1e-4, 
                       fast=as.logical(fast), ada=as.logical(adaptative))
cur.beta.train <- model.refit$beta

### Pour l'evaluation du modele, il faudra calculer :
### 1) la déviance
### 2) le log rank test

##3) PREDICT 
pred <- istacox.predict(model.train=model.refit, x.test=X.test, y.test=y.test)
model <- pred$est

##4) SCORE: get and accumulate the score
pred.score <- istacox.score(y.test=as.matrix(y.test), y.hat=pred$est, type=cv_metric)$perf

save(pred.score, lambda.opt, model, model.refit, file=outputfile)
cat("Writing", outputfile, "\n")