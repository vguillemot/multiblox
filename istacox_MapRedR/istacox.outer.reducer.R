# reducer for outer CV in istacox
args <- commandArgs(trailingOnly = TRUE)

red_file_pattern <- args[1]
nf <- strtoi(args[2])
output_file <- args[3]

global.pred.score <- matrix(NA, nrow=nf, ncol=1)
cur.pred.score <- NULL

for (i in 1:nf) {
  load(paste(red_file_pattern, i,".Rdata", sep=""))
  cur.pred.score[[i]] <- pred.score
}
global.pred.score[, 1] <- as.matrix(cur.pred.score)

print(global.pred.score)

global.stats <- matrix(NA, nrow=3, ncol=1)
global.stats[1,] <- apply(global.pred.score, 2, mean)
global.stats[2,] <- apply(global.pred.score, 2, median)
global.stats[3,] <- apply(global.pred.score, 2, sd)
rownames(global.stats) <- c("Mean","Median","SD")

print(global.pred.score)
print(global.stats)

save(global.pred.score, global.stats, file=output_file)
cat("Writing", output_file, "\n")
