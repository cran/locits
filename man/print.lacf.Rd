\name{print.lacf}
\alias{print.lacf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print lacf class object
}
\description{Prints information about lacf class object.
}
\usage{
\method{print}{lacf}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The lacf class object you want to print
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
}
\author{
Guy Nason
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{lacf}}, \code{\link{plot.lacf}}, \code{\link{summary.lacf}}
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
vlacf <- lacf(v, lag.max=30)
#
# Now let's print the lacf object 
#
print(vlacf)
#Class 'lacf' : Localized Autocovariance/correlation Object:
#       ~~~~  : List with 3 components with names
#             lacf lacr date 
#
#
#summary(.):
#----------
#Name of originating time series:  
#Date produced:  Thu Oct 25 12:11:29 2012 
#Number of times:  256 
#Number of lags:  30 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
