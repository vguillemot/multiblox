# reducer for outter CV in multiblox comparison
args <- commandArgs(trailingOnly = TRUE)

coxnet_red_file_pattern <- args[1]
red_file_pattern <- args[2]
nf <- strtoi(args[3])
output_file <- args[4]

sp <- strsplit(red_file_pattern, "/")
pathtofile <- sp[[1]][1]
for (d in 2:(length(sp[[1]])-1)) {
  pathtofile <- paste(pathtofile, sp[[1]][d], sep="/")
}
pathtofile <- paste(pathtofile, "", sep="/")

load(paste(coxnet_red_file_pattern, 1,".Rdata", sep="")) # on load le 1er fichier juste pour récupérer B
deviance.coxnet <- matrix(NA, nrow=nf, ncol=ncol(global.deviance.coxnet))
pi.coxnet <- matrix(NA, nrow=nf, ncol=ncol(global.pi.coxnet))
coxnet.lambda.opt <- NULL
coxnet.lambda.perblock <- matrix(NA, nrow=nf, ncol=3*B)

designs <- c("hierarchical", "complete")
deviance.multiblox[["hierarchical"]] <- matrix(NA, nrow=nf, ncol=2)
deviance.multiblox[["complete"]] <- matrix(NA, nrow=nf, ncol=2)
pi.multiblox[["hierarchical"]] <- matrix(NA, nrow=nf, ncol=2)
pi.multiblox[["complete"]] <- matrix(NA, nrow=nf, ncol=2)
multiblox.lambda.opt[["hierarchical"]] <- matrix(NA, nrow=nf, ncol=2*B)
multiblox.lambda.opt[["complete"]] <- matrix(NA, nrow=nf, ncol=2*B)
# global.compar <- NULL
deviance.comp <- matrix(NA, nrow=nf, ncol=length(designs))
pi.comp <- matrix(NA, nrow=nf, ncol=length(designs))
labels <- c("LASSO concat", "RIDGE concat", "ElasticNet concat")

for (i in 1:nf) {
  ### deviance.lasso, deviance.ridge, deviance.en, pi.lasso, pi.ridge, pi.en, 
  ### deviance.lasso.perblock, deviance.ridge.perblock, deviance.en.perblock, pi.lasso.perblock, pi.ridge.perblock, pi.en.perblock, 
  ### lambda.opt.lasso, lambda.opt.lasso.perblock, lambda.opt.ridge, lambda.opt.ridge.perblock, 
  ### lambda.opt.en, lambda.opt.en.perblock, B, 
  load(paste(coxnet_red_file_pattern, i,".Rdata", sep="")) ## results for coxnet
  for (j in 1:B) {
      deviance.coxnet[i, ] <- cbind(deviance.coxnet[i, ], deviance.lasso.perblock[1, j], deviance.ridge.perblock[1, j], deviance.en.perblock[1, j])
      pi.coxnet[i, ] <- cbind(pi.coxnet[i, ], pi.lasso.perblock[1, j], pi.ridge.perblock[1, j], pi.en.perblock[1, j])
    }
    deviance.coxnet[i, 3*B+1:3*(B+1)] <- cbind(deviance.lasso, deviance.ridge, deviance.en)
    pi.coxnet[i, 3*B+1:3*(B+1)] <- cbind(pi.lasso, pi.ridge, pi.en)
  }
  coxnet.lambda.opt[[i]] <- cbind(lambda.opt.lasso, lambda.opt.ridge, lambda.opt.en)
  coxnet.lambda.perblock[i, ] <- cbind(lambda.opt.lasso.perblock, lambda.opt.ridge.perblock, lambda.opt.en.perblock)
  ### dev.score, pi.score, lambda.opt, model.refit, design
  for (des in designs){
    load(paste(red_file_pattern, des, i,".Rdata", sep=""))## results for multiblox
    deviance.multiblox[[des]][i, ] <- global.dev.score
    pi.multiblox[[des]][i, ] <- global.pi.score
    multiblox.lambda.opt[[des]][i, ] <- multiblox.lambda.opt 
  }
}
deviance.comp <- cbind(deviance.coxnet, global.dev.score)
pi.comp <- cbind(pi.coxnet, global.pi.score)

for (j in B:1) {
  labels <- c(paste("LASSO block", j , sep=" "), paste("RIDGE block", j , sep=" "), paste("ElasticNet block", j , sep=" "), 
              labels, c("multiblox hierarchical", "multiblox complete")) ## chaque bloc + concat + multiblox
}
print(labels)


# for (cpt in 1:length(designs)) {
#   if (designs[cpt] == "concatenated"){
#     deviance.conf[, B+cpt] <- coxnet.pred.score.concat
#     cat("Optimal lambdas for ", designs[cpt], ":\n     \t", "concat\t\t", sep="")
#     for (i in 1:B) {
#       cat("Block", i, "\t")
#       labels <- c(paste("block", (B - i) + 1, sep=" "), labels)
#     }
#     for (i in 1:nf) {
#       cat(paste("\n[", i, "]\t", sep="") , coxnet.lambda.opt[[i]], "\t")
#       for (j in 1:B) {
#         global.pred.score[i, j] = coxnet.global.pred.score[i, j]
#         cat(coxnet.lambda.perblock[i, j], "\t")
#       }
#     }
#     cat("\n\n", sep="")
#   } else {
#     global.lambda.opt <- matrix(NA, nrow=nf, ncol=2)  # for the plot when B == 2
#     cur.pred.score <- NULL
#     cat("Optimal lambdas for ", designs[cpt], ":\n", sep="")
#     for (i in 1:nf) {
#       load(paste(red_file_pattern, designs[cpt], i,".Rdata", sep=""))
#       cur.pred.score[[i]] <- pred.score
#       global.lambda.opt[i, 1] <- lambda.opt[1]  # for the plot when B == 2
#       global.lambda.opt[i, 2] <- lambda.opt[2]  # for the plot when B == 2
#       cat(paste("[", i, "]\t", sep="") , lambda.opt, "\n", sep=" ")
#     }
#     global.pred.score[, B+cpt] <- as.matrix(cur.pred.score)
#     cat("\n", sep="")
#     # plot global.lambda.opt if B == 2
#     if (B == 2) {
#       pdf(paste(pathtofile, "lambda_", designs[cpt] ,".pdf", sep=""), h=8, w=8)
#       plot(10 ^ (jitter(log10(as.vector(global.lambda.opt[, 1])), amount=0.5)), 
#            10 ^ (jitter(log10(as.vector(global.lambda.opt[, 2])), amount=0.5)), 
#            log="xy", pch="+", xlab="Lambda block 1", ylab="Lambda block2",
#            main="Selected lambdas per fold")
#       dev.off()
#     }
#   }
#   
# }
colnames(deviance.comp) <- colnames(pi.comp) <- labels
# 
# couleurs=list()
# couleurs[["2B"]]=c("yellow", "blue", "green", "red", "blueviolet")
# couleurs[["3B"]]=c("yellow", "blue", "green", "cyan", "red", "blueviolet")
# if(ncol(global.pred.score)==5){
#   coul=couleurs[[1]]  
# } else {
#   coul=couleurs[[2]]
# }
pdf(paste(pathtofile, "Deviance_methods_design_comparison.pdf", sep=""), h=8, w=8)
par(mar=c(9, 5, 3, 2))
boxplot(deviance.comp,  las=2, ylim=c(0, 1), col=palette(rainbow(14)), cex.axis=1.5, cex.lab=1.5, cex.main=1.5,
        main="Deviance for compared methods and designs", xlab="", ylab="Deviance", range=0)
dev.off()

pdf(paste(pathtofile, "PI_methods_design_comparison.pdf", sep=""), h=8, w=8)
par(mar=c(9, 5, 3, 2))
boxplot(pi.comp,  las=2, ylim=c(0, 1), col=palette(rainbow(14)), cex.axis=1.5, cex.lab=1.5, cex.main=1.5,
        main="Prognostic index for compared methods and designs", xlab="", ylab="Prognostic index", range=0)
dev.off()


save(deviance.comp, pi.comp, multiblox.lambda.opt, coxnet.lambda.opt, coxnet.lambda.perblock, file=output_file)
cat("Writing", output_file, "\n")

