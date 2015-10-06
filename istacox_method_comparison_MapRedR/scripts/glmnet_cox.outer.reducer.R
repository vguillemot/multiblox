# reducer for glmnet for Cox model
args <- commandArgs(trailingOnly = TRUE)

red_file_pattern <- args[1]
nf <- strtoi(args[2])
output_file <- args[3]
# labels <- args[4] #GE, CGH ou concat

### deviance.lasso, deviance.ridge, deviance.en, pi.lasso, pi.ridge, pi.en, 
### deviance.lasso.perblock, deviance.ridge.perblock, deviance.en.perblock, pi.lasso.perblock, pi.ridge.perblock, pi.en.perblock, 
### lambda.opt.lasso, lambda.opt.lasso.perblock, lambda.opt.ridge, lambda.opt.ridge.perblock, 
### lambda.opt.en, lambda.opt.en.perblock, B, 
load(paste(red_file_pattern, 1,".Rdata", sep=""))
deviance.coxnet <- matrix(NA, nrow=nf, ncol=3*(B+1))
pi.coxnet <- matrix(NA, nrow=nf, ncol=3*(B+1))
labels <- c("LASSO concat", "RIDGE concat", "ElasticNet concat")
for (i in 1:nf) {
  load(paste(red_file_pattern, i,".Rdata", sep=""))
  for (j in 1:B) {
    global.deviance.coxnet[i, ] <- cbind(global.deviance.coxnet[i, ], deviance.lasso.perblock[1, j], deviance.ridge.perblock[1, j], deviance.en.perblock[1, j])
    global.pi.coxnet[i, ] <- cbind(global.pi.coxnet[i, ], pi.lasso.perblock[1, j], pi.ridge.perblock[1, j], pi.en.perblock[1, j])
  }
  global.deviance.coxnet[i, 3*B+1:3*(B+1)] <- cbind(deviance.lasso, deviance.ridge, deviance.en)
  global.pi.coxnet[i, 3*B+1:3*(B+1)] <- cbind(pi.lasso, pi.ridge, pi.en)
}
for (j in B:1) {
  labels <- c(paste("LASSO block", j , sep=" "), paste("RIDGE block", j , sep=" "), paste("ElasticNet block", j , sep=" "), labels) ## chaque bloc + concat
}
print(labels)

# global.stats <- matrix(NA, nrow=3, ncol=B+1)
# global.stats[1,] <- apply(global.pred.score, 2, mean)
# global.stats[2,] <- apply(global.pred.score, 2, median)
# global.stats[3,] <- apply(global.pred.score, 2, sd)
# rownames(global.stats) <- c("Mean","Median","SD")
# colnames(global.stats) <- labels
colnames(global.deviance.coxnet) <- labels
colnames(global.pi.coxnet) <- labels
print(summary(global.deviance.coxnet))
print(summary(global.pi.coxnet))

pdf("coxnet_performances.pdf")
boxplot(global.deviance.coxnet, main="Deviance for Coxnet models")
boxplot(-log(global.pi.coxnet), main="Prognostic index for Coxnet models", ylab="-log(Pronostic Index)")
dev.off()

save(global.deviance.coxnet, global.pi.coxnet, lambda.opt.lasso, lambda.opt.lasso.perblock, lambda.opt.ridge, lambda.opt.ridge.perblock, 
     lambda.opt.en, lambda.opt.en.perblock, B, file=output_file)
cat("Writing", output_file, "\n")
