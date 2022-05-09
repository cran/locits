\name{summary.lacfCI}
\alias{summary.lacfCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produce a brief summary of the contents of a \code{lacfCI} object
}
\description{Produces brief summary of the contents of a \code{lacfCI}
object.
}
\usage{
\method{summary}{lacfCI}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	The \code{lacfCI} object that you wish to glean info on
}
  \item{\dots}{
	Other arguments.
}
}
\value{No value
}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.
	\doi{10.1111/rssb.12015}}
\author{
Guy Nason.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{print.lacfCI}},\code{\link{Rvarlacf}}
}
\examples{
#
# See example in the Rvarlacf function
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts} 
