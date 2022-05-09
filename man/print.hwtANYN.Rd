\name{print.hwtANYN}
\alias{print.hwtANYN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print out a \code{hwtANYN} class object, eg from the \code{link{hwt}}
	function.
}
\description{Prints out very basic information on an object that
	represents a Haar wavelet transform of a data set of non-dyadic length. 
}
\usage{
\method{print}{hwtANYN}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	The object you wish to print.
}
  \item{\dots}{
	Other arguments
}
}
\value{This function calls the \code{\link{summary.hwtANYN}}
	function as its last action. So, the return from this
	function is the return from \code{\link{summary.hwtANYN}}
}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.
	\doi{10.1111/rssb.12015}}
\author{
Guy Nason.
}

\seealso{\code{\link{hwt}}, \code{\link{summary.hwtANYN}}
}
\examples{
#
# Generate test vector of length 5
#
v2 <- rnorm(5)
#
# Compute Haar wavelet transform
#
v2.hwt <- hwt(v2)
#
# Print out the answer
#
print(v2.hwt)
#Class 'hwtANYN' : Haar Wavelet for Arbitrary Length Data object:
#       ~~~~~~~  : List with 5 components with names
#		   c d nlevels type reindex 
#
#
#summary(.):
#----------
#Levels:  2 
#Filter was: Haar
#Transform type:  wavelet 
#Object was reindex to match wd:  FALSE 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
