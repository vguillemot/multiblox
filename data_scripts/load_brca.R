load.brca.data <-
  function(n=297, pathtofile="./"){
    if (is.null(n)) { n <- 297 }
    #     if (cgh_mode != "norm") {
    #       cgh_mode <- ""
    #     }
    load(paste(pathtofile, "BRCA_TCGA_data_multiblox_297samples.Rdata", sep=""))
    ### X, y
    Xlist <- list(GE=X[["GE"]], MET=X[["MET"]])
    my_var_names <- NULL
    my_var_names[[1]] <- colnames(Xlist[["GE"]])
    my_var_names[[2]] <- colnames(Xlist[["MET"]])
    cat("Block 1 :", nrow(Xlist[[1]]), "x", ncol(Xlist[[1]]), "\n")
    cat("Block 2 :", nrow(Xlist[[2]]), "x", ncol(Xlist[[2]]), "\n")
    cat("Outcome :", dim(y), "\n")
    return (list(X=Xlist, y=as.matrix(y), my_var_names=my_var_names))
  }

read.data <-
  function(cfg){
    pathtodata <- cfg$pathtodata
    n <- cfg$brca$n
    if(length(cfg$brca) == 0) {
      print("bad data configuration")
      quit("no")
    }
    return (load.brca.data(n=n, pathtofile=pathtodata))
  }