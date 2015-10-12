\name{generate_multiblock_survival_simu_data}
\alias{generate_multiblock_survival_simu_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generate_multiblock_survival_simu_data(n = n, p = p, B = B, p_frac = rep(0.25, B), corr = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{p_frac}{
%%     ~~Describe \code{p_frac} here~~
}
  \item{corr}{
%%     ~~Describe \code{corr} here~~
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
function (n = n, p = p, B = B, p_frac = rep(0.25, B), corr = 0.5) 
{
    require(survJamda)
    set.seed(4256)
    ds <- X <- my_var_names <- NULL
    pp <- sum(p)
    ff <- sum(p * p_frac)
    ds <- generate.survival.data(gene.nb = ff, tot.genes = pp, 
        correlation = corr, sample.nb = n, beta.init = 1.5, shape = 1, 
        scale = 1)
    var <- rep(1:B, p * p_frac)
    no_var <- rep(1:B, p - p * p_frac)
    for (b in 1:B) {
        my_var <- ds$ds1[, which(var == b)]
        my_no_var <- ds$ds1[, which(no_var == b)]
        X[[b]] <- cbind(my_var, my_no_var)
        print(dim(X[[b]]))
        print(p[b])
        my_var_names[[b]] <- sprintf("blk\%dvar\%d", b, 1:p[b])
        colnames(X[[b]]) <- my_var_names[[b]]
    }
    return(list(X = X, y = as.matrix(cbind(ds$T, ds$censor)), 
        my_var_names = my_var_names, h0 = 1))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
