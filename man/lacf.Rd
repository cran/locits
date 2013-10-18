\name{lacf}
\alias{lacf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute localized autocovariance.
}
\description{Compute localized autocovariance function for nonstationary
	time series. Note: this function is borrowed from the \code{costat}
	package, and modified to have linear smoothing,  and when that package is complete, it will be removed
	from this package.
}
\usage{
lacf(x, filter.number = 10, family = "DaubLeAsymm", smooth.dev = var,
    AutoReflect = TRUE, lag.max = NULL, WPsmooth.type = "RM",
    binwidth, tol=0.1, maxits=5, ABBverbose=0, verbose=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The time series you wish to analyze
}
  \item{filter.number}{Wavelet filter number you wish to use to 
	analyse the time series (to form the wavelet periodogram, etc)
	See \code{filter.select} for more details.
}
  \item{family}{Wavelet family to use, see \code{filter.select} for
	more details.
}
  \item{smooth.dev}{Change variance estimate for smoothing. Note: \code{var}
	is good for this purpose.
}
  \item{AutoReflect}{If \code{TRUE} then an internal reflection method
	is used to repackage the time series so that it can be
	analyzed by the periodic-assuming wavelet transforms.
}
  \item{lag.max}{The maximum lag of acf required. If NULL then the
	same default as in the regular \code{acf} function is used.
}
  \item{WPsmooth.type}{The type of smoothing used to produce the
	estimate. See \code{\link{ewspec3}} for more advice on this.
}
  \item{binwidth}{If necessary, the \code{binwidth} for the
	spectral smoothing, see \code{\link{ewspec3}} for more info.
	If \code{WTsmooth.type=="RM"} then this argument specifies
	the binwidth of the kernel smoother applied to the wavelet
	periodogram. If the argument is missing or zero then 
	an automatic bandwidth is calculated by \code{\link{AutoBestBW}}.
}
  \item{tol}{Tolerance argument for \code{\link{AutoBestBW}}}
  \item{maxits}{Maximum iterations argument for \code{\link{AutoBestBW}}}
  \item{ABBverbose}{Verbosity of execution of \code{\link{AutoBestBW}}}
  \item{verbose}{If \code{TRUE} then informative message is printed}
  \item{\dots}{Other arguments for \code{\link{ewspec3}}.
}
}
\details{In essence, this routine is fairly simple. First, the EWS
	of the time series is computed. Then formula (14) from
	Nason, von Sachs and Kroisandr (2000) is applied to obtain
	the time-localized autocovariance from the spectral estimate.
}
\value{An object of class \code{lacf} which contains the
	autocovariance. This object can be handled by functions
	from the \code{costat} package. The idea in this package
	is that the function gets used internally and much of the
	same functionality can be achieved by running
	\code{\link{Rvarlacf}} and \code{\link{plot.lacfCI}}. However,
	running \code{lacf} on its own is much faster than
	\code{\link{Rvarlacf}} as the CI computation is intenstive.
}
\references{
Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.

Nason, G.P., von Sachs, R. and Kroisandt, G. (2000) Wavelet processes
	and adaptive estimation of the evolutionary wavelet spectrum.
	\emph{J. R. Statist. Soc.} Ser B, \bold{62}, 271-292.
}
\author{
Guy Nason.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{Rvarlacf}}
}
\examples{
#
# With wavethresh attached, note binwidth is fabricated here,
# just to make the example work. The lacf implementation in
# the costat package performs wavelet (ie maybe better) smoothing automatically
#
v <- lacf(rnorm(256), binwidth=40)
#
# With costat attached also
#
\dontrun{plot(v)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
