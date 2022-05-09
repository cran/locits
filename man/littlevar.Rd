\name{littlevar}
\alias{littlevar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Subsidiary helper function for hwtos2 
}
\description{Computes a variance estimate for \code{\link{hwtos2}}
	Merely takes a wavelet periodogram (actually \code{wd} class
	object), and a level argument. Then extracts the wavelet periodogram
	coefficients at that level and returns twice the mean of their
	squares.
}
\usage{
littlevar(WP, ll)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WP}{The wavelet periodogram that you wish to analyze (actually
	a \code{wd} class object, \code{type="station"}
}
  \item{ll}{A valid level for the periodogram
}
}
\value{Twice the mean of the square of the coefficients at the level
	extracted.
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

\seealso{
\code{\link{hwtos2}}
}
\examples{
#
# Not intended for direct user use
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
