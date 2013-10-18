\name{runmean}
\alias{runmean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute a running mean of a vector
}
\description{This function essentially uses the \code{running.mean}
	function from the \code{igraph} package. However, adjustments
	are made to ensure that the output is always the same length
	as the input (by fiddling at the boundaries).
}
\usage{
runmean(x, binwidth)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector that you wish to smooth using a running mean.
}
  \item{binwidth}{Number of ordinates over which you wish to average
}
}
\details{For example, if \code{binwidth=2} and \code{x=1:6} then
	the function averages each pair to get 1.5, 2.5, 3.5, 4.5, 5.5.
	However, this is only 5 numbers and the input had 6. So, in this
	case the function arranges for the output to be extended (in this
	case 1 gets padded onto the front. For vectors of length > 3
	the padding depends on whether the vector is even or odd.
}
\value{The running mean of the input at the given bandwidth.
}

\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.}

\author{
Guy Nason.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{ewspec3}}
}
\examples{
runmean(1:6, 2)
#
# [1] 1.0 1.5 2.5 3.5 4.5 5.5
#
runmean(1:14, 4)
#
#  [1]  1.75  2.50  3.50  4.50  5.50  6.50  7.50  8.50  9.50 10.50 11.50 12.5
# [13] 13.25 13.50
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
