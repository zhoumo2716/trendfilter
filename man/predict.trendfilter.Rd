% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-trendfilter.R
\name{predict.trendfilter}
\alias{predict.trendfilter}
\title{Predict with trendfilter at new (interior) design points}
\usage{
\method{predict}{trendfilter}(object, newx = NULL, lambda = NULL, deriv = 0L, ...)
}
\arguments{
\item{object}{the result of \code{\link[=trendfilter]{trendfilter()}}}

\item{newx}{numeric vector of new design points at which to evaluate the. The
default, \code{NULL} returns the estimates at the original \code{x} values.}

\item{lambda}{Vector. A user supplied sequence of tuning parameters which
determines the balance between data fidelity and smoothness of the
estimated curve; larger \code{lambda} results in a smoother estimate. The
default, \code{NULL} results in an automatic computation based on \code{nlambda},
the largest value of \code{lambda} that would result in a maximally smooth
estimate, and \code{lambda_min_ratio}. Supplying a value of \code{lambda} overrides
this behaviour. It is likely better to supply a
decreasing sequence of \code{lambda} values than a single (small) value. If
supplied, the user-defined \code{lambda} sequence is automatically sorted in
decreasing order.}

\item{deriv}{integer; the order of the derivative to be evaluated. Default is 0.}

\item{...}{not used}
}
\value{
a vector or matrix with rows corresponding to \code{newx} and columns
corresponding to \code{lambda}
}
\description{
Predict with trendfilter at new (interior) design points
}
\examples{
x <- 1:100 / 101 * 2 * pi
y <- sin(x) + .2 * rnorm(100)
out <- trendfilter(y, x, nlambda = 20L)
predict(out, newx = 1:6, lambda = out$lambda[10])
predict(out, newx = 1:6)
}
