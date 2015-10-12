\name{istacox.predict}
\alias{istacox.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
istacox.predict(model, x, y, D = NULL, lambda, type = c("spll", "deviance", "pi"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
function (model, x, y, D = NULL, lambda, type = c("spll", "deviance", 
    "pi")) 
{
    library(survival)
    source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/sparse.partial.loglik.R")
    source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/partial.loglik.R")
    source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/link.R")
    B <- length(x)
    spll <- res <- lk <- pll <- dev_null <- NULL
    null_model <- matrix(0, nrow = length(model), ncol = 1)
    total_link <- 0
    if (!is.null(D)) {
        for (b in 1:B) {
            lk[[b]] <- link(X = x, D = D, b = b, beta.init = model)
        }
        total_link <- sum(unlist(lk))
    }
    if (type == "spll") {
        print("sparse loglik")
        if (!is.null(D)) {
            for (b in 1:B) {
                spll[[b]] <- sparse.partial.loglik(model = model[[b]], 
                  newdata = x[[b]], newy = y, lambda = lambda[b])
            }
            res <- sum(unlist(spll)) - total_link
        }
        else {
            spll <- sparse.partial.loglik(model = model, newdata = x, 
                newy = y, lambda = lambda)
            res <- spll
        }
    }
    else if (type == "deviance") {
        print("deviance")
        if (!is.null(D)) {
            for (b in 1:B) {
                pll[[b]] <- partial.loglik(model = model[[b]], 
                  newdata = x[[b]], newy = y)
                dev_null[[b]] <- partial.loglik(model = matrix(0, 
                  nrow = length(model[[b]]), ncol = 1), newdata = x[[b]], 
                  newy = y)
            }
            total_pll <- sum(unlist(pll)) - total_link
            total_null <- sum(unlist(dev_null))
            res <- -2 * (total_pll - total_null)
        }
        else {
            pll <- partial.loglik(model = model, newdata = x, 
                newy = y)
            dev_null <- partial.loglik(model = matrix(0, nrow = length(model), 
                ncol = 1), newdata = x, newy = y)
            res <- -2 * (pll - dev_null)
        }
    }
    else if (type == "pi") {
        print("pronostic index")
        if (!is.null(D)) {
            for (b in 1:B) {
                dat[[b]] <- x[[b]] \%*\% model[[b]]
            }
            lp <- Reduce(cbind, x = dat)
            fit <- lapply(coxph(Surv(y[, 1], y[, 2]) ~ lp))
            res <- summary(fit)$logtest["pvalue"]
        }
        else {
            dat <- as.matrix(x) \%*\% as.matrix(model)
            fit <- coxph(Surv(y[, 1], y[, 2]) ~ dat)
            res <- summary(fit)$logtest["pvalue"]
        }
    }
    else {
        print("Sorry, this type of score is not yet implemented !")
    }
    return(list(est = res))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
