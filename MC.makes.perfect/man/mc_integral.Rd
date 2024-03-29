% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{mc_integral}
\alias{mc_integral}
\title{(Improved) Monte Carlo integration}
\usage{
mc_integral(ftn, a, b, n)
}
\arguments{
\item{a}{lower bound of the integration. Must be finite.}

\item{b}{upper bound of the integration. Must be finite and assume b > a.}

\item{n}{number of samples used in the estimation.}
}
\description{
Monte Carlo integral of ftn over [a, b] using a sample of size n.
}
\details{
When we say a Monte Carlo technique is better than another we mean that using the same number of function calls, it has smaller variance.Because our estimates are based on random samples they are themselves random variables.
}
\examples{
f <- function(x) return(x^3 + 1)
mc_integral(f, 0, 1, 10000)
}
\references{
Owen Jones,Robert Maillardet,Andrew Robinson (2014).Introduction to Scientific Programming and Simulation Using R ; Second Edition
}
\seealso{
\code{integral}
}
