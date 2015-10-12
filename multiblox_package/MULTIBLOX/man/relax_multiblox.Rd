\name{relax_multiblox}
\alias{relax_multiblox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
relax_multiblox(x, I, R, D, lambda = 0, eps = 0.001, max.iter = 10000, beta.init = NULL, fast = fast, ada = ada)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{I}{
%%     ~~Describe \code{I} here~~
}
  \item{R}{
%%     ~~Describe \code{R} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{max.iter}{
%%     ~~Describe \code{max.iter} here~~
}
  \item{beta.init}{
%%     ~~Describe \code{beta.init} here~~
}
  \item{fast}{
%%     ~~Describe \code{fast} here~~
}
  \item{ada}{
%%     ~~Describe \code{ada} here~~
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
function (x, I, R, D, lambda = 0, eps = 0.001, max.iter = 10000, 
    beta.init = NULL, fast = fast, ada = ada) 
{
    source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/istacox.R")
    source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/link.R")
    B <- length(x)
    n <- nrow(x[[1]])
    p <- sapply(x, ncol)
    beta <- beta_new <- e <- S <- eta <- cvg <- NULL
    iter <- 1
    cvg[[iter]] <- 10
    beta <- list()
    if (!is.null(beta.init)) {
        beta <- beta_new <- beta.init
    }
    else {
        f <- function(x) x/100
        for (c in 1:B) {
            beta[[c]] <- rep(0, p[c])
        }
        beta_new <- beta
    }
    iter.inner <- 0
    div.inner <- 0
    max.iter.inner <- 200
    consec_max.iter <- 0
    beta_new <- beta
    print(D)
    for (iter in 1:max.iter) {
        for (b in 1:B) {
            print(paste("block : ", b))
            link <- link(x, D, b, beta)
            istacox_res <- istacox(X = x[[b]], I = I, R = R, 
                alpha = 0.5 * lambda, kmax = 1000, epsilon = 1e-04, 
                fast = fast, ada = ada, link = link, beta_init = beta[[b]])
            iter.inner <- iter.inner + istacox_res$k
            print(paste("relax iter : ", iter, sep = ""))
            print(paste("istacox iter : ", iter.inner, sep = ""))
            if (max.iter.inner == istacox_res$k) {
                div.inner <- div.inner + 1
                consec_max.iter <- consec_max.iter + 1
            }
            else {
                consec_max.iter <- 0
            }
            beta_new[[b]] <- istacox_res$beta
        }
        d <- mapply("-", beta, beta_new, SIMPLIFY = FALSE)
        e <- sapply(d, base::norm, "f")
        if (max(e) < eps) 
            break
        beta <- beta_new
        if (consec_max.iter >= (2 * B)) 
            break
    }
    print(paste("Block relaxation iter (inner - div): ", iter, 
        "(", iter.inner, "-", div.inner, ")"))
    return(list(beta = beta_new, convergence = unlist(e), niter = iter))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
