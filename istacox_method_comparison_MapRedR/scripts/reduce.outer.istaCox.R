# reducer for outer CV in istacox
args <- commandArgs(trailingOnly = TRUE)

red_file_pattern <- args[1]
nf <- strtoi(args[2])
output_file <- args[3]

global.dev.score <- matrix(NA, nrow=nf, ncol=1)
global.pi.score <- matrix(NA, nrow=nf, ncol=1)
cur.dev.score <- NULL
cur.pi.score <- NULL

for (i in 1:nf) {
  load(paste(red_file_pattern, i,".Rdata", sep=""))
  cur.dev.score[[i]] <- dev.score
  cur.pi.score[[i]] <- pi.score
}
global.dev.score[, 1] <- as.matrix(cur.dev.score)
global.pi.score[, 1] <- as.matrix(cur.pi.score)

print(global.dev.score)
print(global.pi.score)

global.dev.stats <- matrix(NA, nrow=3, ncol=1)
global.pi.stats <- matrix(NA, nrow=3, ncol=1)

pdf(paste(output_file, ".pdf", sep=""))
boxplot(global.dev.score, main=paste("Deviance for optimal lambda = ", lambda.opt, sep=""), ylab="deviance to null model")
boxplot(-log(global.pi.score), main=paste("Prognostic index for optimal lambda = ", lambda.opt, sep=""), ylab="-log(p.value)")
dev.off()

save(global.dev.score, global.pi.score, file=output_file)
cat("Writing", output_file, "\n")
