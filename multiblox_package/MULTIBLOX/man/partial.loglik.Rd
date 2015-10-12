\name{partial.loglik}
\alias{partial.loglik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
partial.loglik(model, newdata, newy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
}
  \item{newy}{
%%     ~~Describe \code{newy} here~~
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
function (model, newdata, newy) 
{
    beta <- matrix(model, ncol = 1)
    x.o <- newdata[order(newy[, 1]), ]
    y.o <- as.data.frame(newy[order(newy[, 1]), ])
    I <- which(y.o$status == 1)
    R <- lapply(which(y.o$status == 1), function(i) which(y.o$time >= 
        y.o$time[i]))
    names(R) <- paste0("R", which(y.o$status == 1))
    pll <- sum(mapply(function(i) {
        newdata[i, ] \%*\% beta - log(sum(exp(as.matrix(newdata[R[[sprintf("R\%i", 
            i)]], ]) \%*\% beta)))
    }, I))
    return(pll = pll)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
