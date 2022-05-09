\name{hwtos2}
\alias{hwtos2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test of second-order stationarity using wavelets.
}
\description{The main function to perform a test of second-order
	stationarity as outlined in Nason (2012). Essentially,
	this routine computes an evolutionary wavelet spectral
	estimate and then computes the Haar wavelet coefficients
	of each scale of the spectral estimate. Any large Haar
	coefficients are indicative of nonstationarity. A multiple
	hypothesis test assesses whether any of the Haar coefficients
	are large enough to reject the null hypothesis of stationarity.
}
\usage{
hwtos2(x, alpha = 0.05, filter.number = 1, family = "DaubExPhase",
    lowlev = 3, WTscale = NULL, maxSD = NULL, verbose = FALSE,
    silent = FALSE, UseCForVarip2 = TRUE, OPLENGTH = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The time series you want to test for second order
	stationarity. This should be a stochastic series. The function
	will report an error if \code{x} is a constant function. The
	function might not work properly if it contains a significant trend
	or patches of non-stochastic observations.
}
  \item{alpha}{The overall (nominal) size of the test.
}
  \item{filter.number}{The index number of the wavelet used
	to compute the evolutionary spectral estimate with.
}
  \item{family}{The family of wavelet used to compute
	the evolutionary spectral estimate.
}
  \item{lowlev}{Do not compute Haar wavelet coefficients on
	evolutionary wavelet spectra at level lower than \code{lowlev}. 
}
  \item{WTscale}{The theory of the test shows that the Haar wavelet
	coefficients of the raw wavelet periodogram are asymptotically
	normal as long as the scale of the Haar wavelet is
	`coarse' enough. Roughly, speaking \code{WTscale}
	is internally coded to be the log of the square root of T,
	the length of the series (J/2), but you can set another
	value. 
}
  \item{maxSD}{As part of its execution, this function computes an
	evolutionary wavelet spectral estimate from the time series.
	Since the test is based on the assumption of stationarity, the
	EWS is averaged over time. There will be \eqn{J = \log_2 T}
	scale levels and, if \code{maxSD = NULL} then all of the \eqn{J}
	levels get used for later functions, such as computing the
	variance of Haar wavelet coefficients. This argument
	permits you to restrict the number of coarse scales going into
	further calculations (e.g. removes the coarser scales from
	further examination). Mostly, the default will be fine and maximises
	the use of the available information.
}
  \item{verbose}{If \code{TRUE} then informative error messages are printed.
	If \code{FALSE} they are not.
}
  \item{silent}{If \code{TRUE} then no informative messages are printed.
	If \code{FALSE} then a limited amout of informative is printed.
}
  \item{UseCForVarip2}{If \code{TRUE} then fast C code is use to
	compute wavelet coefficients' variance. If \code{FALSE} then
	R code is used wholly throughout, but the execution will be much
	slower.
}
  \item{OPLENGTH}{The \code{PsiJ} and \code{PsiJmat} routines both
	used preallocated storage. This argument can be provided to
	increase the amount of storage. Note, you should not need to change
	this unless the routine as whole stops and tells you to rerun it
	with increased storage.
}
}
\details{This function looks at the Haar wavelet coefficients
	of an evolutionary wavelet spectrum. This is a modification
	of the principle of von Sachs and Neumann (2000) which
	worked with the Haar wavelet coefficients of a local
	Fourier spectrum.

	See also, the \code{stationarity} test which implements
	the Priestley-Subba Rao (1969) test. This function
	is contained in the \code{fractal} package.
}
\value{An object of class \code{tos}, a list containing the following
	components:
\item{nreject}{The number of FDR rejections}
\item{rejpval}{The p-value associated with FDR rejections}
\item{spvals}{A vector of p-values from all of the tests, sorted in
	ascending order.}
\item{sTS}{A vector of sorted test statistics from all of the tests,
	sorted into the same order as \code{spvals}}
\item{AllTS}{A list containing all of the test statistics.
	The first entry contains test statistics corresponding
	to the coarsest scale, the last entry corresponds to the
	finest scale. Each component in the list is either empty
	(because the scale was omitted because it was less than
	\code{lowlev}) or contains a \code{wd} class object. The
	\code{wd} class object contains the test statistics for
	each Haar wavelet coefficient (not the coefficients).
	Hence, the value of the test statistic for any scale/location
	or level of the wavelet periodogram can easily be extracted.
	}
\item{AllPVal}{As \code{AllTS} except the values stored are
	the p-values, not the test statistics.}
\item{alpha}{The nominal size of the overall hypothesis test.}
\item{x}{The original time series that was analyzed}
\item{xSD}{A vector containing \code{J} levels, which is the
	EWS estimate averaged across time.}

}


\references{
Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.
	\doi{10.1111/rssb.12015}

Priestley, M.B. and Subba Rao (1969) A test for non-stationarity of
	time series. \emph{J. R. Statist. Soc. B}, \bold{31}, 140-149.

von Sachs, R. and Neumann, M.H. (2000) A wavelet-based test for
	stationarity. \emph{J. Time Ser. Anal.}, \bold{21},
	597-613.
}
\author{
Guy Nason.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{varip2}}, \code{stationarity}
}
\examples{
#
# First, test a set of iid Gaussians: should be stationary!
#
hwtos2(rnorm(256))
#  8    7    6    5    4    3  
#Class 'tos' : Stationarity Object :
#       ~~~~  : List with 9 components with names
#              nreject rejpval spvals sTS AllTS AllPVal alpha x xSD 
#
#
#summary(.):
#----------
#There are  186  hypothesis tests altogether
#There were  0  FDR rejects
#No p-values were smaller than the FDR val of:   
#Using Bonferroni rejection p-value is  0.0002688172 
#And there would be  0  rejections.
#
#  NOTE: the summary indicates that nothing was rejected: hence stationary!
#
# Second, example. Concatenated Gaussians with different variances
#
hwtos2(c(rnorm(256), rnorm(256,sd=2)))
#  9    8    7    6    5    4    3  
#Class 'tos' : Stationarity Object :
#       ~~~~  : List with 9 components with names
#              nreject rejpval spvals sTS AllTS AllPVal alpha x xSD 
#
#
#summary(.):
#----------
#There are  441  hypothesis tests altogether
#There were  5  FDR rejects
#The rejection p-value was  3.311237e-06 
#Using Bonferroni rejection p-value is  0.0001133787 
#And there would be  5  rejections.
#Listing FDR rejects... (thanks Y&Y!)
#P: 5 HWTlev:  0  indices on next line...[1] 1
#P: 6 HWTlev:  0  indices on next line...[1] 1
#P: 7 HWTlev:  0  indices on next line...[1] 1
#P: 8 HWTlev:  0  indices on next line...[1] 1
#P: 9 HWTlev:  0  indices on next line...[1] 1
#
# NOTE: This time 5 Haar wavelet coefficients got rejected: hence series
# is not stationary.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
