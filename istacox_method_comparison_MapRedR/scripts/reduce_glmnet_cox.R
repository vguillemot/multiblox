# reducer for glmnet for Cox model
args <- commandArgs(trailingOnly = TRUE)

red_file_pattern <- args[1]
nf <- strtoi(args[2])
output_file <- args[3]
# labels <- args[4] #GE, CGH ou concat

load(paste(red_file_pattern, 1,".Rdata", sep=""))
B = ncol(pred.score.perblock)
global.pred.score <- matrix(NA, nrow=nf, ncol=1+B)
labels <- "concat"
for (i in 1:nf) {
  load(paste(red_file_pattern, i,".Rdata", sep=""))
  global.pred.score[i, B+1] <- pred.score.concat  
  for (j in 1:B) {
    global.pred.score[i, j] = pred.score.perblock[1, j]
  }
}
for (j in B:1) {
  labels <- c(paste("block", j , sep=" "), labels) ## chaque bloc + concat
}

global.stats <- matrix(NA, nrow=3, ncol=B+1)
global.stats[1,] <- apply(global.pred.score, 2, mean)
global.stats[2,] <- apply(global.pred.score, 2, median)
global.stats[3,] <- apply(global.pred.score, 2, sd)
rownames(global.stats) <- c("Mean","Median","SD")
colnames(global.stats) <- labels
colnames(global.pred.score) <- labels
print(global.pred.score)
print(global.stats)

save(global.pred.score, global.stats, file=output_file)
cat("Writing", output_file, "\n")
