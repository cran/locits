\name{summary.hwtANYN}
\alias{summary.hwtANYN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize the results of a Haar wavelet transform object computed
	from an arbitrary length vector. 
	
}
\description{This function summarizes the results of a \code{hwtANYN}
	object that contains the results of a Haar wavelet transform that
	had been carried out on an original vector of (potentially) non-dyadic
	length.
}
\usage{
\method{summary}{hwtANYN}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The object that you want a summary of. The object might be
	the results from, e.g., \code{\link{hwt}} function.} 
  \item{...}{Other arguments}
}
\value{None.}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.
	\doi{10.1111/rssb.12015}}
\author{
Guy Nason.
}

\seealso{\code{\link{hwt}}, \code{\link{print.hwtANYN}}
}
\examples{
#
# See help for print.hwtANYN 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
