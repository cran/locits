\name{HwdS}
\alias{HwdS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Compute the non-decimated Haar wavelet transform without
	using periodic boundary conditions.
}
\description{Function uses the \code{filter} function to achieve its aims.
}
\usage{
HwdS(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of dyadic length that you wish to transform.
}
}
\details{The regular \code{wd} function that can compute the non-decimated
	transform uses different kinds of boundary conditions, which can
	result in coefficients being used multiply for consideration in
	a test of stationarity, and distort results. This function
	only computes Haar coefficients on the data it can, without
	wrapround.
}
\value{An object of class \code{wd} which contains the nondecimated
	Haar transform of the input series, \code{x} without periodic
	boundary conditions (nor interval, nor reflection).
}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.}
\author{
Guy Nason.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ewspecHaarNonPer}},
\code{\link{getridofendNA}}
}
\examples{
#
# Apply Haar transform to Gaussian data
#
HwdS(rnorm(32))
#Class 'wd' : Discrete Wavelet Transform Object:
#       ~~  : List with 8 components with names
#              C D nlevels fl.dbase filter type bc date 
#
#$C and $D are LONG coefficient vectors
#
#Created on : Tue Jul 17 15:14:59 2012 
#Type of decomposition:  station 
#
#summary(.):
#----------
#Levels:  5 
#Length of original:  32 
#Filter was:  Haar wavelet 
#Boundary handling:  periodic 
#Transform type:  station 
#Date:  Tue Jul 17 15:14:59 2012 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{algebra}
