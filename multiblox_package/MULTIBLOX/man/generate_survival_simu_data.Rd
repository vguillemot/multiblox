\name{generate_survival_simu_data}
\alias{generate_survival_simu_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generate_survival_simu_data(n, B = 2, B.CGH = 1, rho = rep(0.7, B - 1), p = rep(100, B), n.expl = rep(100, B), noise_corr = rep(0.5, B), alpha = 0.1, seed = 4257, pathtodata = "./")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{B.CGH}{
%%     ~~Describe \code{B.CGH} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{n.expl}{
%%     ~~Describe \code{n.expl} here~~
}
  \item{noise_corr}{
%%     ~~Describe \code{noise_corr} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{pathtodata}{
%%     ~~Describe \code{pathtodata} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n, B = 2, B.CGH = 1, rho = rep(0.7, B - 1), p = rep(100, 
    B), n.expl = rep(100, B), noise_corr = rep(0.5, B), alpha = 0.1, 
    seed = 4257, pathtodata = "./") 
{
    library(mvtnorm)
    library(glmnet)
    if (B.CGH != 0) {
        source(paste(pathtodata, "simulCGHprofiles.R", sep = ""))
    }
    source(paste(pathtodata, "rmvnorm_withGivenLatent.R", sep = ""))
    source(paste(pathtodata, "functions.R", sep = ""))
    set.seed(seed)
    B.norm <- B - B.CGH
    p.CGH <- p[(B.norm + 1):B]
    if (B.CGH != 0) {
        simu.CGH <- NULL
        pc1 <- matrix(0, nrow = n, ncol = B.CGH)
        for (i in 1:B.CGH) {
            simu.CGH[[i]] <- t(as.matrix(simulCGHprofiles(M = n, 
                n = p.CGH[i], k.mean = 8, SNR = 5, lambda = 1, 
                seed = seed)$Y))
            pc1[, i] <- prcomp(simu.CGH[[i]])$x[, 1]
        }
    }
    else {
        pc1 <- matrix(0, nrow = n, ncol = B.norm)
        for (i in 1:B.norm) {
            pc1[, i] <- rnorm(n = n)
        }
    }
    if (length(noise_corr) != length(p)) 
        stop("Attention, chaque bloc doit\n                                            avoir une corr\303\251lation pour le bruit !")
    my_block_names <- paste0("x", 1:B)
    my_var_names <- lapply(1:B, function(b) sprintf("blk\%dvar\%d", 
        b, 1:p[b]))
    k <- 1
    if (length(rho) == 1) {
        sigmaComponents <- toeplitz(c(1, rep(rho, B - 1)))
    }
    else if (length(rho) != B * (B - 1)/2) {
        stop("Attention, la longueur de rho est soit 1 soit (B*(B-1))/2 !")
    }
    else {
        sigmaComponents <- matrix(0, nrow = B, ncol = B)
        for (i in 1:(B - 1)) {
            for (j in (i + 1):B) {
                sigmaComponents[i, j] <- rho[k]
                k <- k + 1
                if (j == B) {
                  break
                }
            }
        }
        I <- idn(sigmaComponents)
        sigmaComponents <- sigmaComponents + t(sigmaComponents) + 
            I
    }
    etaks <- rmvnorm_withGivenLatent(n = n, latent = pc1, sigma = sigmaComponents)
    sigmaz <- mapply(noise_corr, p, FUN = function(nc, q) toeplitz(rep(c(1, 
        nc), c(1, q - 1))), SIMPLIFY = FALSE)
    beta <- NULL
    alphas <- NULL
    z <- x <- eta <- NULL
    for (b in 1:B) {
        alphas[[b]] <- rep(c(alpha, 0), c(n.expl[b], p[b] - n.expl[b]))
        z[[b]] <- rmvnorm(n = n, mean = rep(0, p[[b]]), sigma = sigmaz[[b]])
        x[[b]] <- etaks[, b] \%*\% t(alphas[[b]]) + z[[b]]
    }
    if (B.CGH != 0) {
        for (d in (B.norm + 1):B) {
            x[[d]] <- simu.CGH[[d - B.norm]]
        }
    }
    for (c in 1:B) {
        beta[[c]] <- seq(from = -2 * c, to = 2 * c, length.out = p[c]) + 
            rnorm(p[c], 0, c)
        eta[[c]] <- x[[c]] \%*\% beta[[c]]
    }
    combi <- 0 * rnorm(n, 0, 1)
    for (b in 1:B) {
        print(eta[[b]])
        combi <- combi + eta[[b]]
    }
    ep <- 0.5
    tcens = rbinom(n = n, prob = ep, size = 1)
    h0 <- 0.2
    hx = h0 * exp(combi/n)
    ty = rexp(n, hx)
    y = data.frame(time = ty, status = 1 - tcens, stringsAsFactors = F)
    X = list()
    if (B == 1) {
        X[[1]] = x[[1]]
        names(X) = my_block_names[1:B]
        colnames(X[[1]]) <- my_var_names[[1]]
    }
    else {
        X = x
        names(X) = my_block_names[1:B]
        for (i in 1:B) {
            colnames(X[[i]]) <- my_var_names[[i]]
        }
    }
    return(list(X = X, y = as.matrix(y), my_var_names = my_var_names, 
        h0 = h0))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
