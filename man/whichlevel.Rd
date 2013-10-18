\name{whichlevel}
\alias{whichlevel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Helper routine for \code{mkcoef}
}
\description{Helps \code{\link{mkcoef}} by finding out how many
	more levels are required to compute a set of discrete
	wavelets to a given (other) level.
}
\usage{
whichlevel(J, filter.number = 10, family = "DaubLeAsymm")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{J}{The level that \code{\link{mkcoef}} wants to compute to.
}
  \item{filter.number}{The wavelet number (see \code{wd})
}
  \item{family}{The wavelet family (see \code{wd})
}
}
\details{When computing the discrete wavelets up to a given scale
	we use the inverse wavelet transform to do this. However, to
	generate a wavelet within the range of a wavelet decomposition
	you have to use more scales in the inverse wavelet transform
	than first requested. This is because wavelet coefficients at
	the coarsest scales are associated with wavelets whose support
	is greater than the whole extent of the series. Hence, you 
	have to have a larger wavelet transform, with more levels, insert
	a coefficient mid-level to generate a discrete wavelet whose
	support lies entirely within the extent of the series. This
	function figures out what the extra number of levels should be.
}
\value{Simply returns the required number of levels
}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.}
\author{
Guy Nason.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{mkcoef}}
}
\examples{
whichlevel(6)
# [1] 11
#
# E.g. mkcoef wanted to generate 6 levels of discrete wavelets and
# whichlevel tells it that it needs to generate a wavelet transform
# of at least 11 levels.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
