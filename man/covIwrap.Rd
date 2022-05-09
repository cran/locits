\name{covIwrap}
\alias{covIwrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	A wrapper for the covI function.
}
\description{Computation of the \code{\link{covI}} function is
	intensive. This function permits values of \code{covI}
	to be stored in an object, and then if these values are
	requested again the values can be obtained from a store
	rather than being computed from scratch.
}
\usage{
covIwrap(S, m, n, ll, storewrap, P)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Same argument as for \code{\link{covI}}, a spectral
	estimate (for a stationary series).
}
  \item{m}{The same argument as in \code{\link{covI}}.
}
  \item{n}{The same argument as in \code{\link{covI}}.
}
  \item{ll}{The same argument as in \code{\link{covI}}.
}
  \item{storewrap}{A list. On first call to this function the user should
	supply \code{storewrap=NULL}. This causes the function to
	initialize the storage. On every return from this function
	the storewrap component should be extracted from the list
	and then this storewrap component should be resupplied to
	any future calls to this function. In this way the function
	has access to previously computed values.
}
  \item{P}{Same argument as in \code{\link{covI}}. An autocorrelation
	wavelet computed using the \code{PsiJ} function in \code{wavethresh}.
}
}
\value{A list containing the following components:
\item{ans}{The appropriate covariance}
\item{storewrap}{A list containing information about all previously
	computed covariances. This list should be supplied as the
	\code{storewrap} argument to any future calls of this function,
	so if the same covariance is requested it can be returned from
	\code{storewrap} and not computed again.}
}

\details{Note: \code{covIwrap} could be removed from the function tree
	altogether. I.e. \code{\link{varip2}} could call
	\code{\link{covI}} directly. However,
	\code{covIwrap} considerably improves the efficiency of the algorithm
	as it stores intermediate calculations that can be reused rather
	than being computed repeatedly.
}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.
	\doi{10.1111/rssb.12015}}
\author{
Guy Nason.
}
\seealso{
\code{\link{varip2}}, \code{\link{covI}}
}
\examples{
P1 <- PsiJ(-5, filter.number=1, family="DaubExPhase")
#
# First call to covIwrap
#
ans <- covIwrap(S=c(1/2, 1/4, 1/8, 1/16, 1/32), m=1, n=3, ll=5,
    storewrap=NULL, P=P1)
#
# Make sure you keep the storewrap component.
#
my.storewrap <- ans$storewrap
#
# What is the answer?
#
ans$ans
#[1] 0.8430809
#
# Issue next call to covIwrap: but storewrap argument is now the one we stored.
#
ans <- covIwrap(S=c(1/2, 1/4, 1/8, 1/16, 1/32), m=1, n=3, ll=5,
    storewrap=my.storewrap, P=P1)
#
# This call will reuse the stored value. However, if you change any of the
# arguments then the store won't be used.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{kwd2}% __ONLY ONE__ keyword per line
