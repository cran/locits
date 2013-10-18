\name{covI}
\alias{covI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the covariance between two wavelet periodogram
	ordinates at the same scale, but different time locations.
}
\description{Computes \eqn{cov(I_{\ell, m}, I_{\ell, n})} using the formula
	given in Nason (2012) in Theorem 1. Note: one usually should
	use the \code{\link{covIwrap}} function for efficiency.
}
\usage{
covI(II, m, n, ll, ThePsiJ)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{II}{Actually the *spectral* estimate S, not the periodogram
	values. This is for an assumed stationary series, so this is just
	a vector of length J, one for each scale of S.
}
  \item{m}{Time location m
}
  \item{n}{Time location n
}
  \item{ll}{Scale of the raw wavelet periodogram
}
  \item{ThePsiJ}{Autocorrelation wavelet corresponding to the
	wavelet that computed the raw peridogram (also assumed
	to underlie the time series}
}
\value{The covariance is returned.
}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.}
\author{
Guy Nason.
}

\seealso{\code{\link{covIwrap}}
}
\examples{
P1 <- PsiJ(-5, filter.number=1, family="DaubExPhase")
#
# Compute the covariance
#
covI(II=c(1/2, 1/4, 1/8, 1/16, 1/32), m=1, n=3, ll=5, ThePsiJ=P1)
#
# [1] 0.8430809
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
