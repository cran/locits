\name{AutoBestBW}
\alias{AutoBestBW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Choose a good bandwidth for running mean smoothing of a EWS
spectral estimator.
}
\description{
Computes running mean estimator closest to wavelet estimator of
evolutionary wavelet spectrum.
The idea is to obtain a good linear bandwidth.
}
\usage{
AutoBestBW(x, filter.number = 1, family = "DaubExPhase",
    smooth.dev = var, AutoReflect = TRUE, tol = 0.1, maxits = 5,
    plot.it = FALSE, verbose = 0, ReturnAll = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	Time series you want to analyze.
}
  \item{filter.number}{
	The wavelet filter used to carry out smoothing operations.
}
  \item{family}{
	The wavelet family used to carry out smoothing operations.
}
  \item{smooth.dev}{
	The deviance estimate used for the smoothing (see ewspec help)
}
  \item{AutoReflect}{
	Mitigate periodic boundary conditions of wavelet transforms
	by reflecting time series about RHS end before taking 
	transforms (and is undone before returning the answer).	
}
  \item{tol}{
	Tolerance for golden section search for the best bandwidth
}
  \item{maxits}{
	Maximum number of iterations for the golden section search

}
  \item{plot.it}{
	Plot the values of the bandwidth and its closeness of the
	linear smooth to the wavelet smooth, if TRUE.
	
}
  \item{verbose}{
	If nonzero prints out informative messages about the progress
	of the golden section search. Higher integers produce more
	messages.
}
  \item{ReturnAll}{
	If TRUE then return the best bandwidth (in the ans component),
	the wavelet smooth (in EWS.wavelet) and the closest linear
	smooth (EWS.linear). If FALSE then just the bandwidth is returned.
}
}
\details{Tries to find the best running mean fit to an estimated
	spectrum obtained via wavelet shrinkage. The goal is to try
	and find a reasonable linear bandwidth.
}
\value{If ReturnAll argument is FALSE then the best bandwidth
	is returned.
}

\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.}
\author{
Guy Nason.
}

\seealso{\code{\link{Rvarlacf}}
}
\examples{
#
# Generate synthetic data
#
x <- rnorm(256)
#
# Compute best linear bandwidth
#
tmp <- AutoBestBW(x=x)
#
# Printing it out in my example gives:
# tmp
# [1] 168
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
