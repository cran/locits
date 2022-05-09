\name{print.tosANYN}
\alias{print.tosANYN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print out a \code{tosANYN} class object, eg from the \code{link{hwtos}}
	function.
}
\description{Prints out very basic information on an object that
	represents the output from a test of stationarity.
}
\usage{
\method{print}{tosANYN}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	The object you wish to print.
}
  \item{\dots}{
	Other arguments
}
}
\value{This function calls the \code{\link{summary.tosANYN}}
	function as its last action. So, the return from this
	function is the return from \code{\link{summary.tosANYN}}
}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.
	\doi{10.1111/rssb.12015}}
\author{
Guy Nason.
}

\seealso{\code{\link{hwtos}},\code{\link{summary.tosANYN}}
}
\examples{
#
# See example at end of help for hwtos
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
