\name{link}
\alias{link}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
link(X, D, b, beta.init)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{beta.init}{
%%     ~~Describe \code{beta.init} here~~
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
function (X, D, b, beta.init) 
{
    B <- length(X)
    p <- lapply(X, ncol)
    S <- rep(0, p[[b]])
    if (!is.null(beta.init)) {
        beta <- beta.init
    }
    else {
        f <- function(z) z/100
        for (c in 1:B) {
            beta[[c]] <- rep(0, p[[c]])
        }
    }
    eta <- mapply(X, beta, FUN = function(a, b) a \%*\% b, SIMPLIFY = FALSE)
    for (d in setdiff(1:B, b)) {
        S <- S + D[b, d] * t(X[[b]]) \%*\% eta[[d]]
    }
    return(S)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
