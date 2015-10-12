\name{grad}
\alias{grad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
grad(X, beta, I, R, alpha, link)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{I}{
%%     ~~Describe \code{I} here~~
}
  \item{R}{
%%     ~~Describe \code{R} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{link}{
%%     ~~Describe \code{link} here~~
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
function (X, beta, I, R, alpha, link) 
{
    wij <- mapply(function(i, j) t(exp(X[j, ] \%*\% beta)/sum(exp(X[R[[sprintf("R\%i", 
        i)]], ] \%*\% beta))), I, R)
    names(wij) <- names(R)
    xbar <- t(sapply(names(R), function(r) as.matrix(wij[[r]]) \%*\% 
        as.matrix(X[R[[r]], ])))
    grad <- -colSums(X[I, ] - xbar) + (1 - alpha) * beta - link
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
