\name{fistacox_step_line_search}
\alias{fistacox_step_line_search}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fistacox_step_line_search(X, u, I, R, t = 10, tau = 0.95, alpha, link, kmax = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{u}{
%%     ~~Describe \code{u} here~~
}
  \item{I}{
%%     ~~Describe \code{I} here~~
}
  \item{R}{
%%     ~~Describe \code{R} here~~
}
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{link}{
%%     ~~Describe \code{link} here~~
}
  \item{kmax}{
%%     ~~Describe \code{kmax} here~~
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
function (X, u, I, R, t = 10, tau = 0.95, alpha, link, kmax = 1000) 
{
    Ru_l <- R_ell(X, u, I, R, alpha)
    gradu <- grad(X, u, I, R, alpha, link)
    x <- prox(u - t * gradu, t, alpha)
    Rx_l <- R_ell(X, x, I, R, alpha)
    for (k in 1:kmax) {
        if (Rx_l <= Ru_l + crossprod(gradu, x - u) + 1/(2 * t) * 
            sum((x - u)^2)) 
            break
        t <- tau * t
        x <- prox(u - t * gradu, t, alpha)
        Rx_l <- R_ell(X, x, I, R, alpha)
    }
    return(t)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
