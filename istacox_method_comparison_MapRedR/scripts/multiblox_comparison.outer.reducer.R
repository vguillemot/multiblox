# reducer for outter CV in multiblox comparison
args <- commandArgs(trailingOnly = TRUE)

glm_red_file_pattern <- args[1]
red_file_pattern <- args[2]
nf <- strtoi(args[3])
output_file <- args[4]

sp <- strsplit(red_file_pattern, "/")
pathtofile <- sp[[1]][1]
for (d in 2:(length(sp[[1]])-1)) {
  pathtofile <- paste(pathtofile, sp[[1]][d], sep="/")
}
pathtofile <- paste(pathtofile, "", sep="/")

load(paste(glm_red_file_pattern, 1,".Rdata", sep="")) # on load le 1er fichier juste pour récupérer B
B = ncol(pred.score.perblock)
glm.global.pred.score <- matrix(NA, nrow=nf, ncol=B)
labels <- "concat"
glm.pred.score.concat <- rep(NA, nf)
glm.lambda.opt <- NULL
glm.lambda.perblock <- matrix(NA, nrow=nf, ncol=B)
for (i in 1:nf) {
  load(paste(glm_red_file_pattern, i,".Rdata", sep=""))
  glm.pred.score.concat[i] <- pred.score.concat 
  glm.lambda.opt[[i]] <- lambda.opt
  glm.lambda.perblock[i, ] <- lambda.opt.perblock
  for (j in 1:B) {
    glm.global.pred.score[i, j] = pred.score.perblock[1, j]
  }
}
for (j in B:1) {
  labels <- c(paste("block", j , sep=" "), labels) ## chaque bloc + concat
}

designs <- c("concatenated", "hierarchical", "complete")
labels <- designs
# global.compar <- NULL
global.pred.score <- matrix(NA, nrow=nf, ncol=length(designs)+B)

for (cpt in 1:length(designs)) {
  if (designs[cpt] == "concatenated"){
    global.pred.score[, B+cpt] <- glm.pred.score.concat
    cat("Optimal lambdas for ", designs[cpt], ":\n     \t", "concat\t\t", sep="")
    for (i in 1:B) {
      cat("Block", i, "\t")
      labels <- c(paste("block", (B - i) + 1, sep=" "), labels)
    }
    for (i in 1:nf) {
      cat(paste("\n[", i, "]\t", sep="") , glm.lambda.opt[[i]], "\t")
      for (j in 1:B) {
        global.pred.score[i, j] = glm.global.pred.score[i, j]
        cat(glm.lambda.perblock[i, j], "\t")
      }
    }
    cat("\n\n", sep="")
  } else {
    global.lambda.opt <- matrix(NA, nrow=nf, ncol=2)  # for the plot when B == 2
    cur.pred.score <- NULL
    cat("Optimal lambdas for ", designs[cpt], ":\n", sep="")
    for (i in 1:nf) {
      load(paste(red_file_pattern, designs[cpt], i,".Rdata", sep=""))
      cur.pred.score[[i]] <- pred.score
      global.lambda.opt[i, 1] <- lambda.opt[1]  # for the plot when B == 2
      global.lambda.opt[i, 2] <- lambda.opt[2]  # for the plot when B == 2
      cat(paste("[", i, "]\t", sep="") , lambda.opt, "\n", sep=" ")
    }
    global.pred.score[, B+cpt] <- as.matrix(cur.pred.score)
    cat("\n", sep="")
    # plot global.lambda.opt if B == 2
    if (B == 2) {
      pdf(paste(pathtofile, "lambda_", designs[cpt] ,".pdf", sep=""), h=8, w=8)
      plot(10 ^ (jitter(log10(as.vector(global.lambda.opt[, 1])), amount=0.5)), 
           10 ^ (jitter(log10(as.vector(global.lambda.opt[, 2])), amount=0.5)), 
           log="xy", pch="+", xlab="Lambda block 1", ylab="lambda block2",
           main="lambdas per fold")
      dev.off()
    }
  }
  
}
colnames(global.pred.score) <- labels

global.stats <- matrix(NA, nrow=3, ncol=length(labels))
global.stats[1,] <- apply(global.pred.score, 2, mean)
global.stats[2,] <- apply(global.pred.score, 2, median)
global.stats[3,] <- apply(global.pred.score, 2, sd)
colnames(global.stats) <- labels
rownames(global.stats) <- c("Mean","Median","SD")

print(global.pred.score)
print(global.stats)
couleurs=list()
couleurs[["2B"]]=c("yellow", "blue", "green", "red", "blueviolet")
couleurs[["3B"]]=c("yellow", "blue", "green", "cyan", "red", "blueviolet")
if(ncol(global.pred.score)==5){
  coul=couleurs[[1]]  
} else {
  coul=couleurs[[2]]
}
pdf(paste(pathtofile, "exactConcordanceIndex_methods_design_comparison.pdf", sep=""), h=8, w=8)
par(mar=c(9, 5, 3, 2))
boxplot(global.pred.score,  las=2, ylim=c(0, 1), col=coul, cex.axis=1.5, cex.lab=1.5, cex.main=1.5,
        main="Exact Concordance Index for compared methods and designs", xlab="", ylab="Exact Concordance Index", range=0)
dev.off()


save(global.pred.score, global.stats, file=output_file)
cat("Writing", output_file, "\n")

