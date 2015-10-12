\name{istacox_step_line_search}
\alias{istacox_step_line_search}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
istacox_step_line_search(X, beta, I, R, t = 10, tau = 0.5, alpha, link, kmax = 1000)
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
function (X, beta, I, R, t = 10, tau = 0.5, alpha, link, kmax = 1000) 
{
    R_l <- R_ell(X, beta, I, R, alpha)
    g <- grad(X, beta, I, R, alpha, link)
    Gt <- neg_step(X, t, beta, alpha, g)
    Rt_l <- R_ell(X, beta - t * Gt, I, R, alpha)
    for (k in 1:kmax) {
        if (Rt_l <= R_l - t * t(g) \%*\% Gt + t/2 * sum(Gt^2)) 
            break
        t <- tau * t
        Gt <- neg_step(X, t, beta, alpha, g)
        Rt_l <- R_ell(X, beta - t * Gt, I, R, alpha)
    }
    return(t)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
