# Mapper for glmnet Cox
library(glmnet)
library(survival)
#library(multiblog)
source("functions.R")
source("multiblox.score.R")

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
method <- cv_method
outer_cv_metric <- cv_metric
glm.model.concat <- res.glm.concat <- NULL
pred.score.concat <- matrix(NA, nrow=1, ncol=1)
lambda.opt <- matrix(NA, nrow=1, ncol=1)
lglmnet_cv.perblock <- NULL
pred.score.perblock <- matrix(NA, nrow=1, ncol=B)
lambda.opt.perblock <- matrix(NA, nrow=1, ncol=B)

# scaling
if (scale) {
  scl_fun <- function(data, scaled) {
    scale(data, center = attr(scaled, "scaled:center"),
          scale = attr(scaled, "scaled:scale")) }
  X.train.concat <- scale2(X.train.concat)
  X.test.concat <- scl_fun(data=X.test.concat, scaled=X.train.concat)
}

# n.zeros.test <- length(y.test[which(y.test==0)])
# n.uns.test <- length(y.test[which(y.test==1)])

# concat
lglmnet_cv <- cv.glmnet(x = X.train.concat, y = Surv(y.train[, 1], y.train[, 2]), family=family, alpha = 0, 
                             type.measure=type.measure)
attach(what=as.data.frame(X.train.concat))
y.pred <- predict(lglmnet_cv, newx=X.test.concat, s="lambda.min")
detach(name=as.data.frame(X.train.concat))
pred.score.concat[1] <- multiblox.score(as.numeric(y.pred), as.matrix(y.test), metric=outer_cv_metric)$perf
lambda.opt[1] <- lglmnet_cv$lambda.min
# per block
lglmnet_cv.perblock <- y.chapo.perblock <- NULL
for (j in 1:B) {
  lglmnet_cv.perblock[[j]] <- cv.glmnet(x = X.train[[j]], y = y.train, family=family, alpha = 0, 
                                          type.measure=type.measure)
  attach(what=as.data.frame(X.train[[j]]))
  y.pred <- predict(lglmnet_cv.perblock[[j]], newx=X.test[[j]], s="lambda.min")
  detach(name=as.data.frame(X.train[[j]]))
  pred.score.perblock[1, j] <- multiblox.score(as.numeric(y.pred), as.matrix(y.test), metric=outer_cv_metric)$perf
  lambda.opt.perblock[1, j] <- lglmnet_cv.perblock[[j]]$lambda.min
}

print(pred.score.concat)
print(pred.score.perblock)

save(pred.score.concat, pred.score.perblock, lambda.opt, lambda.opt.perblock, B, file=outputfile)
cat("Writing", outputfile, "\n")
