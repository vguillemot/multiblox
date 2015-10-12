\name{read.data}
\alias{read.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.data(cfg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cfg}{
%%     ~~Describe \code{cfg} here~~
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
function (cfg) 
{
    pathtodata <- cfg$pathtodata
    B <- cfg$simu$B
    B.CGH <- cfg$simu$B.CGH
    seed <- cfg$simu$seed
    rho <- cfg$simu$rho
    n <- cfg$simu$n
    p <- cfg$simu$p
    n.expl <- cfg$simu$n.expl
    alpha <- cfg$simu$alpha
    noise_corr <- cfg$simu$noise_corr
    if (length(cfg$simu) == 0) {
        print("bad data configuration")
        quit("no")
    }
    if ((is.null(B)) + (is.null(rho)) + (is.null(seed)) + (is.null(n)) + 
        (is.null(alpha)) + (is.null(p)) + (is.null(noise_corr)) + 
        (length(p) != B) + (length(noise_corr) != B)) {
        print("Missing or bad parameters to generate simulated data")
        quit("no")
    }
    res <- generate_survival_simu_data(n = n, B = B, B.CGH = B.CGH, 
        rho = rho, p = p, n.expl = n.expl, noise_corr = noise_corr, 
        alpha = alpha, seed = seed, pathtodata = pathtodata)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
