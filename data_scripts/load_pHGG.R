load.pHGG.data <-
  function(cgh_mode=c("seg", "norm"), n=92, pathtofile="./"){
    if (is.null(n)) { n <- 92 }
    if (is.null(cgh_mode)) { cgh_mode <- "seg" }
#     if (cgh_mode != "norm") {
#       cgh_mode <- ""
#     }
    load(paste(pathtofile, "pHGG_multiblox_data.Rdata", sep=""))
### X, y, clinic, CGH_annot
    my_var_names <- NULL
    my_var_names[[1]] <- colnames(X[["GE"]])
    my_var_names[[2]] <- colnames(X[["CGH"]])
    cat("Block 1 :", nrow(X[[1]]), "x", ncol(X[[1]]), "\n")
    cat("Block 2 :", nrow(X[[2]]), "x", ncol(X[[2]]), "\n")
    cat("Outcome :", dim(y), "\n")
    return (list(X=X, y=as.matrix(y), my_var_names=my_var_names))
  }

read.data <-
  function(cfg){
    pathtodata <- cfg$pathtodata
    n <- cfg$pHGG$n
    cgh_mode <- cfg$pHGG$cgh_mode
    if(length(cfg$pHGG) == 0) {
      print("bad data configuration")
      quit("no")
    }
    return (load.pHGG.data(cgh_mode=cgh_mode, n=n, pathtofile=pathtodata))
  }