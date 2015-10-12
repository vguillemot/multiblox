\name{multiblox.lambda.tune}
\alias{multiblox.lambda.tune}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
multiblox.lambda.tune(X, y, D, trainmat, i, outer_it, lambda.grid, scale = T, method = "CV", metric = "spll", adaptative = TRUE, fast = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
}
  \item{trainmat}{
%%     ~~Describe \code{trainmat} here~~
}
  \item{i}{
%%     ~~Describe \code{i} here~~
}
  \item{outer_it}{
%%     ~~Describe \code{outer_it} here~~
}
  \item{lambda.grid}{
%%     ~~Describe \code{lambda.grid} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{metric}{
%%     ~~Describe \code{metric} here~~
}
  \item{adaptative}{
%%     ~~Describe \code{adaptative} here~~
}
  \item{fast}{
%%     ~~Describe \code{fast} here~~
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
function (X, y, D, trainmat, i, outer_it, lambda.grid, scale = T, 
    method = "CV", metric = "spll", adaptative = TRUE, fast = TRUE) 
{
    cat("### outer CV : ", outer_it, ", inner CV iteration : ", 
        i, "\n")
    library(MASS)
    library(CMA)
    library(survival)
    source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/relax_multiblox.R")
    source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/istacox.predict.R")
    source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/istacox.score.R")
    source("/home/philippe/github/multiblox/istacox_method_comparison_MapRedR/scripts/functions.R")
    B <- length(X)
    if (B == 1) {
        N <- nrow(X[[1]])
    }
    else {
        N <- nrow(X)
    }
    print(D)
    beta.train <- eta.train <- eta.test <- NULL
    pred.score <- res <- model <- cox.model <- cv <- NULL
    beta0 <- NULL
    for (b in 1:B) {
        beta0[[b]] <- matrix(0, nrow = (ncol(X[[b]])), ncol = 1)
    }
    ind <- trainmat[i, ]
    y.train <- y[ind, ]
    y.test <- y[-ind, ]
    X.train <- lapply(X, function(mm) mm[ind, ])
    if (method == "LOOCV") {
        X.test <- lapply(X, function(mm) t(as.matrix(mm[-ind, 
            ])))
    }
    else {
        X.test <- lapply(X, function(mm) {
            as.matrix(mm[-ind, ])
        })
    }
    if (scale) {
        X.train <- lapply(X.train, function(mm) scale2(mm))
        scl_fun <- function(data, scaled) {
            scale(data, center = attr(scaled, "scaled:center"), 
                scale = attr(scaled, "scaled:scale"))
        }
        X.test <- mapply(scl_fun, X.test, X.train, SIMPLIFY = FALSE)
    }
    x.o <- list()
    x.o <- lapply(X.train, function(b) b[order(y.train[, 1]), 
        ])
    y.o <- as.data.frame(y.train[order(y.train[, 1]), ])
    colnames(y.o) <- c("time", "status")
    I.train <- which(y.o$status == 1)
    R.train <- lapply(which(y.o$status == 1), function(i) which(y.o$time >= 
        y.o$time[i]))
    names(R.train) <- paste0("R", which(y.o$status == 1))
    for (l in 1:nrow(lambda.grid)) {
        cat("lambdas : ", l, " -> ")
        for (bid in 1:length(lambda.grid[l, ])) {
            cat(" ", lambda.grid[l, bid], sep = "")
        }
        cat("\n")
        print(adaptative)
        res[[l]] <- relax_multiblox(x = x.o, I.train, R.train, 
            D = D, lambda = lambda.grid[l, ], max.iter = 1000, 
            ada = as.logical(adaptative), fast = as.logical(fast), 
            beta.init = beta0)
        beta.train[[l]] <- res[[l]]$beta
        pred <- istacox.predict(model = beta.train[[l]], D = D, 
            x = x.o, y = y.o, lambda = lambda.grid[l, ], type = metric)
        pred.score[[l]] <- istacox.score(as.matrix(y.test), pred[["est"]])[["perf"]]
        beta0 <- res[[l]][["beta"]]
    }
    return(list(res = res, pred.score = pred.score))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
