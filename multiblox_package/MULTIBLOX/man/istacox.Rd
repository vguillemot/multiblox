\name{istacox}
\alias{istacox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
istacox(X, I, R, alpha, kmax = 1000, epsilon = 1e-04, fast = FALSE, ada = FALSE, link, beta_init)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
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
  \item{kmax}{
%%     ~~Describe \code{kmax} here~~
}
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
}
  \item{fast}{
%%     ~~Describe \code{fast} here~~
}
  \item{ada}{
%%     ~~Describe \code{ada} here~~
}
  \item{link}{
%%     ~~Describe \code{link} here~~
}
  \item{beta_init}{
%%     ~~Describe \code{beta_init} here~~
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
function (X, I, R, alpha, kmax = 1000, epsilon = 1e-04, fast = FALSE, 
    ada = FALSE, link, beta_init) 
{
    p <- ncol(X)
    n <- nrow(X)
    if (!is.null(beta_init)) {
        betaold <- beta_init
    }
    else {
        betaold <- rnorm(p)
    }
    t <- 1/max(eigen(t(X) \%*\% X)$values)
    betanew <- prox(betaold - t * grad(X, betaold, I, R, alpha, 
        link), t, alpha)
    for (k in 2:kmax) {
        if (fast) {
            u <- betanew + (k - 1)/(k + 2) * (betanew - betaold)
        }
        else {
            u <- betanew
        }
        if (ada) {
            if (fast) {
                t <- fistacox_step_line_search(X, u, I, R, t, 
                  tau = 0.95, alpha, link)
            }
            else {
                t <- istacox_step_line_search(X, betaold, I, 
                  R, t, tau = 0.95, alpha, link)
            }
        }
        betaold <- betanew
        betanew <- prox(u - t * grad(X, u, I, R, alpha, link), 
            t, alpha)
        if (sum((betanew - betaold)^2) < epsilon) 
            break
    }
    return(list(beta = betanew, k = k))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
