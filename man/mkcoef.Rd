\name{mkcoef}
\alias{mkcoef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute discrete wavelets.
}
\description{For a given wavelet computes a list with each entry
	of the list containing that discrete wavelet at a different
	scale. The first entry corresponds to the finest wavelet,
	the next entry to the next finest, and so on.
}
\usage{
mkcoef(J, filter.number = 10, family = "DaubLeAsymm")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{J}{A NEGATIVE integer. -J is the maximum number of levels to compute.
}
  \item{filter.number}{The filter number (number of vanishing moments)
	of the underlying wavelet to use.
}
  \item{family}{The family of the wavelet. See \code{wd} help for further
	info.
}
}
\value{A list of length J. The first entry contains the discrete
	wavelet at the finest scale, the 2nd entry contains the next
	most finest wavelet, and so on.
}

\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.}

\author{
Guy Nason.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{Rvarlacf}},
\code{\link{whichlevel}}
}
\examples{
#
# E.g. compute discrete Haar wavelets on scales 1, 2, 3.
#
mkcoef(-3, 1, "DaubExPhase")
#[[1]]
#[1]  0.7071068 -0.7071068
#
#[[2]]
#[1]  0.5  0.5 -0.5 -0.5
#
#[[3]]
#[1]  0.3535534  0.3535534  0.3535534  0.3535534 -0.3535534 -0.3535534 -0.3535534
#[8] -0.3535534
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{algebra}
