\name{make.lambda.grid}
\alias{make.lambda.grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.lambda.grid(X, path = c("naive", "smart", "norm1"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{path}{
%%     ~~Describe \code{path} here~~
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
function (X, path = c("naive", "smart", "norm1")) 
{
    B <- length(X)
    scale <- TRUE
    if (scale) {
        X <- lapply(X, function(mm) scale2(mm))
    }
    lambda.max <- rep(0, B)
    eps <- 0.005
    for (i in 1:B) {
        lambda.max[[i]] <- 0.5 * sum(apply(X[[i]], 2, function(c) norm(c, 
            "2")))
    }
    lambda.min <- eps * lambda.max
    nb.lambda <- 10 - 1 - 1
    l <- NULL
    j <- c(0, 0:nb.lambda)
    for (i in 1:B) {
        l[[i]] <- lambda.max[[i]] * (lambda.min[[i]]/lambda.max[[i]])^(j/nb.lambda)
        l[[i]][[1]] <- 10000
    }
    if (B == 1) {
        print("un seul bloc !!")
        lg <- as.matrix(expand.grid(l))
        lambda.grid <- matrix(0, nrow = length(lg), ncol = 1)
        lambda.grid[, 1] <- lg
    }
    else {
        lambda.grid <- as.matrix(expand.grid(l))
    }
    if (path == "naive") {
        n1 <- nrow(lambda.grid):1
    }
    else if (path == "norm1") {
        n1 = apply(lambda.grid, 1, function(u) u[[1]] + u[[2]])
    }
    else if (path == "smart") {
        n1 <- nrow(lambda.grid):1
    }
    return(lambda.grid[order(n1, decreasing = T), ])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
