\name{tvar1sim}
\alias{tvar1sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Simulate a realization from a particular TVAR(1) model.
}
\description{Simulates a realization from a TVAR(1) model where
	the AR(1) parameter moves from 0.9 to -0.9 in equal steps
	over 512 time points. The realization is also of length 512.
	The innovations are normally distributed with mean zero and
	standard deviation of \code{sd}. 
}
\usage{
tvar1sim(sd = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sd}{This is the standard deviation of the Gaussian innovation.
}
}
\details{This function is easily converted into one that does the
	same thing but for a different sample size.
}
\value{A realization of the aforementioned TVAR(1) process.
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

\seealso{\code{\link{Rvarlacf}}
}
\examples{
#
# Generate realization from the TVAR(1) process
#
x <- tvar1sim()
#
# Maybe plot it
#
\dontrun{ts.plot(x)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
