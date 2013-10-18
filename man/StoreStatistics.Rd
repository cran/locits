\name{StoreStatistics}
\alias{StoreStatistics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Interogates calculation store to see how well we are
	reusing previous calculations (debugging)
}
\description{The computation of the variance of the lacf
	estimator is intensive and we try to speed that up
	by reusing calculations. These calculations are stored
	in internal C arrays. This function interrogates those
	arrays and can provide details on how well the storage
	is working and provide hints if more storage needs to be
	allocated. For very large time series it is possible that
	values need to be calculated that can be stored and this
	function can monitor this.
}
\usage{
StoreStatistics()
}
%- maybe also 'usage' for other objects documented here.
\details{The function prints out the state of the storage.
	Three numbers are reported on. 1. The number of values
	that were calculated but not stored "outside framework".
	Ideally you want this number to be low, if it gets persistently
	high then more storage needs to be allocated in the
	C code (notably MAXELL, MAXJ, MAXK, MAXD for the
	ThmStore and ValExists arrays). 

	The other two numbers are "Number stored" and "Number found".
	The first number corresponds to the number of values calculated
	once and then stored. The second number contains the
	number of times the software interogated the store and found
	a value that it did not have to then calculate. So,
	ideally, you'd like the latter number to be a high percentage
	of the former number, as this means the store is working
	efficiently.

	Note, this function is definitely not intended for
	casual users. However, for users of very large series,
	who have the computational resources, these storage parameters
	might need to be increased.

	The values will be zero if \code{\link{Rvarlacf}} has not
	yet been called, and only refer to the last call to that function
	(as the function zeroes the store on invocation).
}
\value{None.
}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.}
\author{
Guy Nason.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{Rvarlacf}}
}
\examples{
#
# Simulate some data
#
x <- tvar1sim()
#
# Calculate lacf and confidence intervals
#
x.lacf <- Rvarlacf(x=x, nz=50, var.lag.max=20)
#
# Find out how the store did
#
StoreStatistics()
#Number calculated outside framework:  0 
#Number calculated then stored:  154440 
#Number found in store:  14980680 
#
#Overall % calculated:  1.020408 
#% outside framework:  0 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
