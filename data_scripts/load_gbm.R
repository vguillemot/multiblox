load.gbm.data <-
  function(n=369, pathtofile="./"){
    if (is.null(n)) { n <- 369 }
    #     if (cgh_mode != "norm") {
    #       cgh_mode <- ""
    #     }
    load(paste(pathtofile, "GBM_TCGA_data_multiblox_369samples.Rdata", sep=""))
    ### X, y
    Xlist <- list(GE=t(X[["GE"]]), CNASNP=t(X[["CNASNP"]]))
    my_var_names <- NULL
    my_var_names[[1]] <- colnames(Xlist[["GE"]])
    my_var_names[[2]] <- colnames(Xlist[["CNASNP"]])
    cat("Block 1 :", nrow(Xlist[[1]]), "x", ncol(Xlist[[1]]), "\n")
    cat("Block 2 :", nrow(Xlist[[2]]), "x", ncol(Xlist[[2]]), "\n")
    cat("Outcome :", dim(y), "\n")
    return (list(X=Xlist, y=as.matrix(y), my_var_names=my_var_names))
  }

read.data <-
  function(cfg){
    pathtodata <- cfg$pathtodata
    n <- cfg$gbm$n
    if(length(cfg$gbm) == 0) {
      print("bad data configuration")
      quit("no")
    }
    return (load.gbm.data(n=n, pathtofile=pathtodata))
  }