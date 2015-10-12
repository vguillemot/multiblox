\name{simulCGHprofiles}
\alias{simulCGHprofiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simulCGHprofiles(M, n, k.mean, SNR, lambda, seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{k.mean}{
%%     ~~Describe \code{k.mean} here~~
}
  \item{SNR}{
%%     ~~Describe \code{SNR} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
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
function (M, n, k.mean, SNR, lambda, seed) 
{
    set.seed(seed)
    P = c(0.5, 0, 0.5)
    Pcum = cumsum(P)
    K = 1 + rpois(M, k.mean - 1)
    k.max = max(K)
    X = matrix(nrow = M, ncol = n)
    E = matrix(nrow = M, ncol = n)
    mu = matrix(nrow = M, ncol = n)
    tau = matrix(NaN, nrow = M, ncol = k.max)
    for (m in (1:M)) {
        tau.tmp = ceiling((n - 1) * runif(K[m] - 1))
        tau.tmp = sort(tau.tmp)
        tau[m, (1:K[m])] = c(tau.tmp, n)
        level = vector(length = K[m])
        level[1] = 0
        mu.tmp = rep(level[1], tau[m, 1])
        if (K[m] > 1) {
            for (k in (2:K[m])) {
                if (mu.tmp[tau[m, k - 1]] != 0) 
                  level[k] = 0
                else level[k] = -1 + sum(runif(1) > Pcum)
                mu.tmp = c(mu.tmp, rep(level[k], (tau[m, k] - 
                  tau[m, k - 1])))
            }
        }
        mu[m, ] = mu.tmp
    }
    mu[mu == 1] = 0.58
    sigma0 = sqrt(sum(mu^2)/(n * M * SNR))
    TrueRef = t(matrix(rnorm(n, 0, sigma0/lambda) + (sigma0/lambda) * 
        sin(2 * pi * (1:n)/n * 10), ncol = M, nrow = n))
    for (m in (1:M)) {
        E[m, ] = rnorm(n, 0, sigma0)
    }
    X = mu + TrueRef + E
    Y = data.frame(t(X))
    colnames(Y) = paste(rep("Ind", M), c(1:M), sep = "")
    invisible(list(Y = Y, Km = K, mu0 = mu, sigma0 = sigma0, 
        theta0 = TrueRef[1, ]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
