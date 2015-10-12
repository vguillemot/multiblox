\name{rmvnorm_withGivenLatent}
\alias{rmvnorm_withGivenLatent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rmvnorm_withGivenLatent(n, latent, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{latent}{
%%     ~~Describe \code{latent} here~~
}
  \item{mean}{
%%     ~~Describe \code{mean} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{pre0.9_9994}{
%%     ~~Describe \code{pre0.9_9994} here~~
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
function (n, latent, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
    method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE) 
{
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
        check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    sigma1 <- sigma
    dimnames(sigma1) <- NULL
    if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
        warning("sigma is numerically not symmetric")
    }
    method <- match.arg(method)
    if (method == "eigen") {
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- ev$vectors \%*\% diag(sqrt(ev$values), length(ev$values)) \%*\% 
            t(ev$vectors)
    }
    else if (method == "svd") {
        sigsvd <- svd(sigma)
        if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- t(sigsvd$v \%*\% (t(sigsvd$u) * sqrt(sigsvd$d)))
    }
    else if (method == "chol") {
        retval <- chol(sigma, pivot = TRUE)
        o <- order(attr(retval, "pivot"))
        retval <- retval[, o]
    }
    lat <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994)
    lat[, (ncol(sigma) - ncol(latent)):ncol(sigma)] <- latent
    retval <- lat \%*\% retval
    retval <- sweep(retval, 2, mean, "+")
    colnames(retval) <- names(mean)
    retval
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
