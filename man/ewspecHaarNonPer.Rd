\name{ewspecHaarNonPer}
\alias{ewspecHaarNonPer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute evolutionary wavelet spectrum (EWS) estimate based
	on the Haar wavelet transform.
}
\description{This function uses the special \code{HwdS} function
	to compute the Haar wavelet transform with out boundary
	conditions (neither periodic, interval, mirror reflection).
	This is so all coefficients are genuine Haar coefficients without
	involving extra/repeated data.

}
\usage{
ewspecHaarNonPer(x, filter.number = 1, family = "DaubExPhase",
    UseLocalSpec = TRUE, DoSWT = TRUE, WPsmooth = TRUE,
    verbose = FALSE, smooth.filter.number = 10,
    smooth.family = "DaubLeAsymm",
    smooth.levels = 3:WPwst$nlevels - 1, smooth.dev = madmad,
    smooth.policy = "LSuniversal", smooth.value = 0,
    smooth.by.level = FALSE, smooth.type = "soft",
    smooth.verbose = FALSE, smooth.cvtol = 0.01,
    smooth.cvnorm = l2norm, smooth.transform = I,
    smooth.inverse = I)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of dyadic length that contains the time series
	you want to form the EWS of.
}
  \item{filter.number}{Should always be 1 (for Haar)
}
  \item{family}{Should always be "DaubExPhase", for Haar.
}
  \item{UseLocalSpec}{
	Should always be \code{TRUE}.
}
  \item{DoSWT}{Should always be \code{TRUE}
}
  \item{WPsmooth}{Should alway be \code{TRUE} to do smoothing. If \code{FALSE}
	then not smoothed.
%%     ~~Describe \code{WPsmooth} here~~
}
  \item{verbose}{ If \code{TRUE} informative messages are printed during
	the progress of the algorithm.
}
  \item{smooth.filter.number}{Wavelet filter number for doing the
	wavelet smoothing of the EWS estimate.
}
  \item{smooth.family}{Wavelet family for doing the wavelet smoothing
	of the EWS estimate.
}
  \item{smooth.levels}{Which levels of the EWS estimate to apply smoothing
	to.
}
  \item{smooth.dev}{What kind of deviance to use. The default is madmad,
	an alternative might be var. 
}
  \item{smooth.policy}{What kind of smoothing to use. See help
	page for \code{ewspec}
}
  \item{smooth.value}{If a manual value has to be supplied according
	to the \code{smooth.policy} then this is it.
}
  \item{smooth.by.level}{If \code{TRUE} then all levels are smoothed
	independently with different smoothing, otherwise all levels
	are smoothed together (eg one threshold for all levels).
}
  \item{smooth.type}{The type of wavelet smoothing "hard" or "soft"
}
  \item{smooth.verbose}{If \code{TRUE} then informative messages about
	the smoothing are printed.
}
  \item{smooth.cvtol}{If cross-validation smoothing is used, this is
	the tolerance
}
  \item{smooth.cvnorm}{If cross-validation smoothing used, this is the
	norm that's used
}
  \item{smooth.transform}{
	A transform is applied before smoothing
}
  \item{smooth.inverse}{
	The inverse transform is applied after smoothing
}
}
\details{
This function is very similar
	to \code{ewspec} from wavethresh, and many arguments here perform
	the same function as there.
}
\value{The same value as for the \code{ewspec} function.
}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.}
\author{
Guy Nason.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{hwtos2}},
\code{\link{HwdS}}
}
\examples{
#
# Requires wavethresh, so not run directly in installation of package
#
ewspecHaarNonPer(rnorm(512))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
