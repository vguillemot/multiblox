\name{generate_monoblock_survival_data}
\alias{generate_monoblock_survival_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generate_monoblock_survival_data(n, p, p_frac)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{p_frac}{
%%     ~~Describe \code{p_frac} here~~
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
function (n, p, p_frac) 
{
    set.seed(10101)
    N <- n
    nzc <- p * p_frac
    x <- matrix(rnorm(N * p), nrow = N, ncol = p)
    betastar <- rnorm(nzc)
    fx <- x[, seq(nzc)] \%*\% (betastar * p_frac)
    h0 <- 1
    hx <- h0 * exp(fx)
    ty <- rexp(N, hx)
    tcens <- rbinom(n = N, prob = 0.3, size = 1)
    y <- cbind(time = ty, status = 1 - tcens)
    X = list()
    X[[1]] = x
    my_var_names <- sprintf("blk1var\%d", 1:p)
    colnames(X[[1]]) <- my_var_names
    return(list(X = X, y = as.matrix(y), my_var_names = my_var_names, 
        h0 = h0))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
