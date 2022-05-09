\name{summary.lacf}
\alias{summary.lacf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizes a lacf object
}
\description{Summarizes a lacf object
}
\usage{
\method{summary}{lacf}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The lacf object you wish summarized.
}
  \item{\dots}{Other arguments
}
}
\value{None
}
\references{
Cardinali, A. and Nason, G.P. (2012) Costationarity of Locally 
Stationary Time Series using \code{costat}.

Cardinali, A. and Nason, G.P. (2010) Costationarity of locally stationary
time series. \emph{J. Time Series Econometrics}, \bold{2}, Issue 2, Article 1.

Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.
	\doi{10.1111/rssb.12015}
}
\author{
Guy Nason
}


\seealso{\code{\link{lacf}}, \code{\link{plot.lacf}}, \code{\link{print.lacf}}
}
\examples{
#
# Make some dummy data, e.g. white noise
#
v <- rnorm(256)
#
# Compute the localized autocovariance (ok, the input is stationary
# but this is just an example. More interesting things could be achieved
# by putting the results of simulating from a LSW process, or piecewise
# stationary by concatenating different stationary realizations, etc.
#
vlacf <- lacf(v, lag.max=20)
#
# Now let's summarize the lacf object
#
summary(vlacf)
#Name of originating time series:  
#Date produced:  Thu Oct 25 12:11:29 2012 
#Number of times:  256 
#Number of lags:  20  
}
\keyword{ts}
