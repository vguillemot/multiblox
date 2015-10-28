# Mapper for Coxnet Ridge, LASSO, ElasticNet
library("glmnet")
library("methods") ### to avoid the following bug in glmnet when executed via Rscript :
# Error in is(x, "CsparseMatrix") : could not find function "new"
# Calls: source ... <Anonymous> -> glmnet -> elnet -> getcoef -> drop0 -> is
# Execution halted
library("survival")
library("MULTIBLOX")
# source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/functions.R")
# source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/istacox.predict.R")

args <- commandArgs(trailingOnly = TRUE)
# common parameters
fold_file <- args[1]
outputfile <- args[2]
# outer_fold <- strtoi(args[3])
# custom parameters
# None

# type <- "class"
type.measure <- "deviance" # deviance of partial likelihood for Cox model
family <- "cox"

### Initialization
load(fold_file)
B <- length(X.train)
colnames(y.train) <- colnames(y.test) <- c("time", "status")
method <- cv_method
outer_cv_metric <- cv_metric

lglmnet_lasso<- lglmnet_ridge <- lglmnet_en <- NULL
deviance.lasso <- deviance.ridge <- deviance.en <- matrix(NA, nrow=1, ncol=1)
pi.lasso <- pi.ridge <- pi.en <- matrix(NA, nrow=1, ncol=1)
lambda.opt.lasso <- lambda.opt.ridge <- lambda.opt.en <- matrix(NA, nrow=1, ncol=1)

lglmnet_lasso.perblock <- lglmnet_ridge.perblock <- lglmnet_en.perblock <- NULL
deviance.lasso.perblock <- deviance.ridge.perblock <- deviance.en.perblock <- matrix(NA, nrow=1, ncol=B)
pi.lasso.perblock <- pi.ridge.perblock <- pi.en.perblock <- matrix(NA, nrow=1, ncol=B)
lambda.opt.lasso.perblock <- lambda.opt.ridge.perblock <- lambda.opt.en.perblock <- matrix(NA, nrow=1, ncol=B)

# scaling
if (scale) {
  scl_fun <- function(data, scaled) {
    scale(data, center = attr(scaled, "scaled:center"),
          scale = attr(scaled, "scaled:scale")) }
  X.train.concat <- scale2(X.train.concat)
  X.test.concat <- scl_fun(data=X.test.concat, scaled=X.train.concat)
}

# concat
print("concat")
### alpha = 0 -> LASSO
lglmnet_lasso <- cv.glmnet(x = X.train.concat, y = Surv(y.train[, 1], y.train[, 2]), family=family, alpha = 0, 
                             type.measure=type.measure, grouped=TRUE)
### alpha = 1 -> Ridge
lglmnet_ridge <- cv.glmnet(x = X.train.concat, y = Surv(y.train[, 1], y.train[, 2]), family=family, alpha = 1, 
                           type.measure=type.measure, grouped=TRUE)
### alpha = NULL -> ElasticNet
lglmnet_en <- cv.glmnet(x = X.train.concat, y = Surv(y.train[, 1], y.train[, 2]), family=family, 
                           type.measure=type.measure, grouped=TRUE)
# attach(what=as.data.frame(X.train.concat))
# y.pred.lasso <- predict(lglmnet_lasso, newx=X.test.concat, s="lambda.min")
# y.pred.ridge <- predict(lglmnet_ridge, newx=X.test.concat, s="lambda.min")
# y.pred.en <- predict(lglmnet_en, newx=X.test.concat, s="lambda.min")
# detach(name=as.data.frame(X.train.concat))

print("deviance")
### deviance
lambda.opt.lasso[1] <- lglmnet_lasso$lambda.min
deviance.lasso[1] <- istacox.predict(model=lglmnet_lasso$glmnet.fit$beta[, which(lglmnet_lasso$lambda==lglmnet_lasso$lambda.min)], x=X.test.concat, y=as.matrix(y.test), D=NULL, lambda=lambda.opt.lasso[1], type="deviance")$est
lambda.opt.ridge[1] <- lglmnet_ridge$lambda.min
deviance.ridge[1] <- istacox.predict(model=lglmnet_ridge$glmnet.fit$beta[, which(lglmnet_ridge$lambda==lglmnet_ridge$lambda.min)], x=X.test.concat, y=as.matrix(y.test), D=NULL, lambda=lambda.opt.ridge[1], type="deviance")$est
lambda.opt.en[1] <- lglmnet_en$lambda.min
deviance.en[1] <- istacox.predict(model=lglmnet_en$glmnet.fit$beta[, which(lglmnet_en$lambda==lglmnet_en$lambda.min)], x=X.test.concat, y=as.matrix(y.test), D=NULL, lambda=lambda.opt.en[1], type="deviance")$est

print("pi")
### pronostic score
pi.lasso[1] <- istacox.predict(model=lglmnet_lasso$glmnet.fit$beta[, which(lglmnet_lasso$lambda==lglmnet_lasso$lambda.min)], x=X.test.concat, y=as.matrix(y.test), D=NULL, lambda=lambda.opt.lasso[1], type="pi")$est
pi.ridge[1] <- istacox.predict(model=lglmnet_ridge$glmnet.fit$beta[, which(lglmnet_ridge$lambda==lglmnet_ridge$lambda.min)], x=X.test.concat, y=as.matrix(y.test), D=NULL, lambda=lambda.opt.ridge[1], type="pi")$est
pi.en[1] <- istacox.predict(model=lglmnet_en$glmnet.fit$beta[, which(lglmnet_en$lambda==lglmnet_en$lambda.min)], x=X.test.concat, y=as.matrix(y.test), D=NULL, lambda=lambda.opt.en[1], type="pi")$est

print("per block")
# per block
lglmnet_lasso.perblock <- lglmnet_ridge.perblock <- lglmnet_en.perblock <- NULL
for (j in 1:B) {
  lglmnet_lasso.perblock[[j]] <- cv.glmnet(x = X.train[[j]], y = y.train, family=family, alpha = 0, 
                                          type.measure=type.measure, grouped=TRUE)
  lglmnet_ridge.perblock[[j]] <- cv.glmnet(x = X.train[[j]], y = y.train, family=family, alpha = 1, 
                                           type.measure=type.measure, grouped=TRUE)
  lglmnet_en.perblock[[j]] <- cv.glmnet(x = X.train[[j]], y = y.train, family=family,
                                           type.measure=type.measure, grouped=TRUE)
#   attach(what=as.data.frame(X.train[[j]]))
#   y.pred.lasso <- predict(lglmnet_lasso.perblock[[j]], newx=X.test[[j]], s="lambda.min")
#   y.pred.ridge <- predict(lglmnet_ridge.perblock[[j]], newx=X.test[[j]], s="lambda.min")
#   y.pred.en <- predict(lglmnet_en.perblock[[j]], newx=X.test[[j]], s="lambda.min")
  #detach(name=as.data.frame(X.train[[j]]))
  lambda.opt.lasso.perblock[1, j] <- lglmnet_lasso.perblock[[j]]$lambda.min
  lambda.opt.ridge.perblock[1, j] <- lglmnet_ridge.perblock[[j]]$lambda.min
  lambda.opt.en.perblock[1, j] <- lglmnet_en.perblock[[j]]$lambda.min
  
  deviance.lasso.perblock[1, j] <- istacox.predict(model=lglmnet_lasso.perblock[[j]]$glmnet.fit$beta[, which(lglmnet_lasso.perblock[[j]]$lambda==lglmnet_lasso.perblock[[j]]$lambda.min)], x=X.test[[j]], y=as.matrix(y.test), D=NULL, lambda=lambda.opt.lasso.perblock[1, j], type="deviance")$est
  deviance.ridge.perblock[1, j] <- istacox.predict(model=lglmnet_ridge.perblock[[j]]$glmnet.fit$beta[, which(lglmnet_ridge.perblock[[j]]$lambda==lglmnet_ridge.perblock[[j]]$lambda.min)], x=X.test[[j]], y=as.matrix(y.test), D=NULL, lambda=lambda.opt.ridge.perblock[1, j], type="deviance")$est
  deviance.en.perblock[1, j] <- istacox.predict(model=lglmnet_en.perblock[[j]]$glmnet.fit$beta[, which(lglmnet_en.perblock[[j]]$lambda==lglmnet_en.perblock[[j]]$lambda.min)], x=X.test[[j]], y=as.matrix(y.test), D=NULL, lambda=lambda.opt.en.perblock[1, j], type="deviance")$est
  pi.lasso.perblock[1, j] <- istacox.predict(model=lglmnet_lasso.perblock[[j]]$glmnet.fit$beta[, which(lglmnet_lasso.perblock[[j]]$lambda==lglmnet_lasso.perblock[[j]]$lambda.min)], x=X.test[[j]], y=as.matrix(y.test), D=NULL, lambda=lambda.opt.lasso.perblock[1, j], type="pi")$est
  pi.ridge.perblock[1, j] <- istacox.predict(model=lglmnet_ridge.perblock[[j]]$glmnet.fit$beta[, which(lglmnet_ridge.perblock[[j]]$lambda==lglmnet_ridge.perblock[[j]]$lambda.min)], x=X.test[[j]], y=as.matrix(y.test), D=NULL, lambda=lambda.opt.ridge.perblock[1, j], type="pi")$est
  pi.en.perblock[1, j] <- istacox.predict(model=lglmnet_en.perblock[[j]]$glmnet.fit$beta[, which(lglmnet_en.perblock[[j]]$lambda==lglmnet_en.perblock[[j]]$lambda.min)], x=X.test[[j]], y=as.matrix(y.test), D=NULL, lambda=lambda.opt.en.perblock[1, j], type="pi")$est
  
}

print(paste("deviance LASSO concat : ", deviance.lasso, sep=""))
print(paste("deviance LASSO per block  : ", deviance.lasso.perblock, sep=""))
print(paste("deviance RIDGE concat : ", deviance.ridge, sep=""))
print(paste("deviance RIDGE per block  : ", deviance.ridge.perblock, sep=""))
print(paste("deviance EN concat : ", deviance.en, sep=""))
print(paste("deviance EN per block  : ", deviance.en.perblock, sep=""))

print(paste("Prognostic Index LASSO concat : ", pi.lasso, sep=""))
print(paste("Prognostic Index LASSO per block  : ", pi.lasso.perblock, sep=""))
print(paste("Prognostic Index RIDGE concat : ", pi.ridge, sep=""))
print(paste("Prognostic Index RIDGE per block  : ", pi.ridge.perblock, sep=""))
print(paste("Prognostic Index EN concat : ", pi.en, sep=""))
print(paste("Prognostic Index EN per block  : ", pi.en.perblock, sep=""))

save(deviance.lasso, deviance.ridge, deviance.en, pi.lasso, pi.ridge, pi.en, 
     deviance.lasso.perblock, deviance.ridge.perblock, deviance.en.perblock, pi.lasso.perblock, pi.ridge.perblock, pi.en.perblock, 
     lambda.opt.lasso, lambda.opt.lasso.perblock, lambda.opt.ridge, lambda.opt.ridge.perblock, 
     lambda.opt.en, lambda.opt.en.perblock, B, file=outputfile)
cat("Writing", outputfile, "\n")
