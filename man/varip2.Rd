\name{varip2}
\alias{varip2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Direct computation of estimate of variance of v_ip,
	the Haar wavelet coefficients of the periodogram.
}
\description{Performs a direct computation of an estimate
	of the variance of the Haar wavelet coefficients of
	the raw wavelet periodogram of a time series. 
}
\usage{
varip2(i, p, ll, S, P)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{Scale parameter of Haar wavelet analyzing periodogram.
	Scale 1 is the finest scale.
}
  \item{p}{Location parameter of Haar wavelet analyzing periodogram
}
  \item{ll}{Scale of the raw wavelet periodogram being analyzed
}
  \item{S}{
	Estimate of the spectrum, under the assumption of stationarity.
	So, this is just a vector of (possibly) J scales (which is often
	the usual spectral estimate averaged over time). Note: that the
	main calling function, \code{\link{hwtos2}}, actually passes
	\code{maxD} levels.
}
  \item{P}{Is an autocorrelation wavelet object, returned by the
	\code{PsiJ} function. The wavelet concerned is the analyzing
	one underlying the raw wavelet periodogram of the series.
}
}
\details{Computes the variance of the Haar wavelet coefficients of
	the raw wavelet periodogram. Note, that this is merely
	an estimate of the variances.
}
\value{A list with the following components:
\item{covAA}{A component of the variance}
\item{covAB}{A component of the variance}
\item{covBB}{A component of the variance}
\item{ans}{The actual variance}
}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.
	\doi{10.1111/rssb.12015}}
\author{
Guy Nason.
}

\seealso{\code{\link{Cvarip2}},\code{\link{hwtos2}}, \code{\link{covIwrap}}
}
\examples{
#
# Generate autocorrelation wavelets
#
P1 <- PsiJ(-5, filter.number=1, family="DaubExPhase")
#
#
# Now compute varip2: this is the variance of the Haar wavelet coefficient
# at the finest scale, location 10 and P1 autocorrelation wavelet.
# Note, I've used S to be the exact coefficients which would be for white noise.
# In practice, S would be an *estimate* calculated from the data.
#
varip2(i=1, p=10, ll=2, S=c(1/2, 1/4, 1/8, 1/16, 1/32), P=P1)
#
# Ans component is 1.865244
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
