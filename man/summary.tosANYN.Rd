\name{summary.tosANYN}
\alias{summary.tosANYN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize the results of a test of stationarity contained
in an \code{tosANYN} class object. 
}
\description{This function summarizes the results of a \code{tosANYN}
	object that contains information about a test of stationarity.
	The summary function prints out information about how many
	individual hypothesis tests there were, how many were rejected.
	If the hypothesis of stationarity
	is rejected then the routine also prints out a list of the Haar
	wavelet coefficients (and their scales, locations and scale of
	the wavelet periodogram) that were significant.
	The function also returns a lot of this information (invisibly).
}
\usage{
\method{summary}{tosANYN}(object, quiet = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The output from a function that carries out a test
	of stationarity, e.g. \code{\link{hwtos}}.}
  \item{quiet}{If this argument is \code{TRUE} then nothing is printed
	to the screen, although the information is still returned as
	an object. If \code{FALSE} then this function prints out
	the information about the hypothesis test.
}
  \item{...}{Other arguments}
}
\value{The function returns a list which contain a list of the
	rejected coefficients. Each list item contains the index
	of a particular rejected coefficient, which is a vector
	of at least three elements. The first element corresponds to
	the scale of the wavelet periodogram, the second is the level of
	the Haar wavelet transform, and all remaining values are the
	index of the significant wavelet coefficients at that Haar wavelet
	transform scale. The list also contains the total number of
	Haar wavelet coefficients rejected and the \code{mctype}
	argument also.
}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.}
\author{
Guy Nason.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{hwtos}}, \code{\link{print.tosANYN}}
}
\examples{
#
# See example for hwtos, this contains two examples where
# summary.tosANYN (as summary(.) is used in the output.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
