load.ov.data <-
  function(n=555, pathtofile="./"){
    if (is.null(n)) { n <- 555 }
    #     if (cgh_mode != "norm") {
    #       cgh_mode <- ""
    #     }
    load(paste(pathtofile, "OV_TCGA_data_multiblox_555samples.Rdata", sep=""))
    ### X, y
    Xlist <- list(GE=t(X[["GE"]]), MET=t(X[["MET"]]))
    my_var_names <- NULL
    my_var_names[[1]] <- colnames(Xlist[["GE"]])
    my_var_names[[2]] <- colnames(Xlist[["MET"]])
    cat("Block 1 :", nrow(Xlist[[1]]), "x", ncol(Xlist[[1]]), "\n")
    cat("Block 2 :", nrow(Xlist[[2]]), "x", ncol(Xlist[[2]]), "\n")
    cat("Outcome :", dim(y), "\n")
    return (list(X=Xlist, y=y, my_var_names=my_var_names))
  }

read.data <-
  function(cfg){
    pathtodata <- cfg$pathtodata
    n <- cfg$ov$n
    if(length(cfg$ov) == 0) {
      print("bad data configuration")
      quit("no")
    }
    return (load.ov.data(n=n, pathtofile=pathtodata))
  }