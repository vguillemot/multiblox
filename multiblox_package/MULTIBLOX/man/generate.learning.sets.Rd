\name{generate.learning.sets}
\alias{generate.learning.sets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generate.learning.sets(N, y, method = "MCCV", nf.cv = 10, inner_cv_seed = 4257)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{nf.cv}{
%%     ~~Describe \code{nf.cv} here~~
}
  \item{inner_cv_seed}{
%%     ~~Describe \code{inner_cv_seed} here~~
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
function (N, y, method = "MCCV", nf.cv = 10, inner_cv_seed = 4257) 
{
    library(CMA)
    set.seed(inner_cv_seed)
    if (nf.cv == N) {
        trainmat <- GenerateLearningsets(N, y, method = "LOOCV")@learnmatrix
        attr(trainmat, "method") <- "LOOCV"
    }
    else {
        if (method == "CV") {
            trainmat <- GenerateLearningsets(N, y, method = "CV", 
                fold = nf.cv, strat = TRUE)@learnmatrix
            attr(trainmat, "method") <- "CV"
        }
        else {
            trainmat <- GenerateLearningsets(N, y, method = "MCCV", 
                niter = nf.cv, ntrain = floor(N * 0.8), strat = TRUE)@learnmatrix
            attr(trainmat, "method") <- "MCCV"
        }
    }
    return(trainmat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
