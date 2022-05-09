\name{summary.tos}
\alias{summary.tos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize the results of a test of stationarity contained i
	a \code{tos} object.
	
}
\description{This function summarizes the results of a \code{tos}
	object that contains information about a test of stationarity.
	The summary function prints out information about how many
	individual hypothesis tests there were, how many were rejected
	and what the (equivalent) rejection p-value was (in the last cases
	both for FDR and Bonferroni). If the hypothesis of stationarity
	is rejected then the routine also prints out a list of the Haar
	wavelet coefficients (and their scales, locations and scale of
	the wavelet periodogram) that were significant.
	The function also returns a lot of this information (invisibly).
}
\usage{
\method{summary}{tos}(object, size = 0.05, quiet = FALSE, mctype = "FDR", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The output from a function that carries out a test
	of stationarity, e.g. \code{\link{hwtos2}}.
}
  \item{size}{The nominal overall significance level of the test.
}
  \item{quiet}{If this argument is \code{TRUE} then nothing is printed
	to the screen, although the information is still returned as
	an object. If \code{FALSE} then this function prints out
	the information about the hypothesis test.
}
  \item{mctype}{This argument can be \code{"FDR"} for false discovery
	rate or \code{"BON"} for Bonferroni. This argument does not
	effect the basis printout. However, it does control whether
	FDR or Bonferroni rejects are listed, and it does control
	the type of information returned by the function (whether
	FDR or Bonferroni info).
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
	\bold{75}, 879-904.
	\doi{10.1111/rssb.12015}}
\author{
Guy Nason.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{hwtos2}}, \code{\link{print.tos}}
}
\examples{
#
# See example for hwtos2, this contains two examples where
# summary.tos (as summary(.) is used in the output.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
