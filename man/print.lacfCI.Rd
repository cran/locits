\name{print.lacfCI}
\alias{print.lacfCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print basic information about a \code{lacfCI} object.
}
\description{Prints basic information about a \code{lacfCI} object, which
	contains information on confidence intervals for localized
	autocovariance.
}
\usage{
\method{print}{lacfCI}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The \code{lacfCI} object.
}
  \item{\dots}{
	Other arguments
}
}
\value{The last action of this function is to compute
	\code{\link{summary.tos}} so the return code is whatever
	that function returns.
}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.}
\author{
Guy Nason.
}

\seealso{\code{\link{summary.lacfCI}}, \code{\link{Rvarlacf}}
}
\examples{
#
# See example on Rvarlacf help page
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts} 
