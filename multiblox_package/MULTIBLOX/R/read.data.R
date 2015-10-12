read.data <-
function(cfg){
    pathtodata <- cfg$pathtodata
    B <- cfg$simu$B
    B.CGH <- cfg$simu$B.CGH
    seed <- cfg$simu$seed
    rho <- cfg$simu$rho
    n <- cfg$simu$n
    p <- cfg$simu$p
    n.expl <- cfg$simu$n.expl
    alpha <- cfg$simu$alpha
    noise_corr <- cfg$simu$noise_corr
    if(length(cfg$simu) == 0) {
      print("bad data configuration")
      quit("no")
    }
    # basic configuration checking
    if ((is.null(B)) + (is.null(rho)) + (is.null(seed)) + (is.null(n)) + (is.null(alpha)) +
          (is.null(p)) + (is.null(noise_corr)) + (length(p) != B) + (length(noise_corr) != B)) {
      print("Missing or bad parameters to generate simulated data")
      quit("no")
    }
    res <- generate_survival_simu_data(n=n, B=B, B.CGH=B.CGH, rho=rho, p=p, n.expl=n.expl, 
                              noise_corr=noise_corr, alpha=alpha, seed=seed, 
                              pathtodata=pathtodata)
    return (res)
  }
