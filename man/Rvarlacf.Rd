\name{Rvarlacf}
\alias{Rvarlacf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Compute confidence intervals for localized autocovariance for
        locally stationary time series.
}
\description{Compute a localized autocovariance and associated
	confidence intervals for a locally stationary time
	series. The underlying theory assumes a locally stationary
	wavelet time series, but will work well for other time
	series that are not too far away.

}
\usage{
Rvarlacf(x, nz, filter.number = 1, family = "DaubExPhase",
    smooth.dev = var, AutoReflect = TRUE, lag.max = NULL,
    WPsmooth.type = "RM", binwidth = 0, mkcoefOBJ, ThePsiJ,
    Cverbose = 0, verbose = 0, OPLENGTH = 10^5, var.lag.max = 3,
    ABB.tol = 0.1, ABB.plot.it = FALSE, ABB.verbose = 0,
    ABB.maxits = 10, truedenom=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The time series you wish to analyze
}
  \item{nz}{The time point at which you wish to compute the
	localized autocovariance for.
}
  \item{filter.number}{The analysis wavelet for many things, including
	smoothing. See \code{wd} for information on the various types.
}
  \item{family}{The analysis wavelet family. See \code{wd} again.
}
  \item{smooth.dev}{The deviance function used to perform smoothing of
	the evolutionary wavelet spectrum. 
}
  \item{AutoReflect}{The internal wavelet transforms assume periodic
	boundary conditions. However, most time series are not
	periodic (in terms of their support, e.g. the series at time
	1 is not normally anywhere near the value of the series at
	time T). This argument, if \code{TRUE} mitigates this
	by reflecting the whole series by the right-hand end, computing
	the transform (for which periodic transforms are now valid)
	and then junking the second half of the estimate. Although this
	is slightly more computationally intensive, the results are better.
}
  \item{lag.max}{The maximum number of lags to compute the localized
	autocovariance for. The default is the same as in the
	regular \code{acf} function.
}
  \item{WPsmooth.type}{The type of smoothing of the evolutionary
	wavelet spectrum and the localized autocovariance. See the
	arguments to \code{\link{lacf}}.
}
  \item{binwidth}{The smoothing bandwidth associated with the
	smoothing controlled by \code{WPsmooth.type}. If this value
	is zero then the \code{binwidth} is computed automatically
	by the routine. And if \code{verbose>0} the value is also
	printed.
}
  \item{mkcoefOBJ}{Optionally, the appropriate discrete wavelet
	transform object can be supplied. If it is not supplied
	then the routine automatically computes it. There is a small
	saving in providing it, so for everyday use probably not worth
	it.  
}
  \item{ThePsiJ}{As for \code{mkcoefOBJ} but the autocorrelation
	wavelet object.
}
  \item{Cverbose}{If positive integer then the called C code produces
	verbose messages. Useful for debugging. 
}
  \item{verbose}{If positive integer >0 then useful messages are printed.
	Higher values give more information.
}
  \item{OPLENGTH}{Parameter that controls storage allocated to
	the \code{PsiJ} routine. It is possible, for large time series,
	you might be asked to increase this value.
}
  \item{var.lag.max}{Number of lags that you want to compute confidence
	intervals for. Usually, it is quick to compute for more lags,
	so this could usually be set to be the value of \code{lag.max}
	above. 
}
  \item{ABB.tol}{
	The routine selects the automatic bandwidth via a golden section
	search. This argument controls the optimization tolerance.
}
  \item{ABB.plot.it}{Whether or not to plot the iterations of the
	automatic bandwidth golden section search. (\code{TRUE}/\code{FALSE})
}
  \item{ABB.verbose}{Positive integer controlling the amount of
	detail from the automatic bandwidth golden section search
	algorithm. If zero nothing is produced.
}
  \item{ABB.maxits}{
	The maximum number of iterations in the automatic bandwidth
	golden section search.
}
  \item{truedenom}{If TRUE use the actual number of terms in the sum as
	the denominator in the formula for the calculation of the
	covariance of the smoothed periodogram. If FALSE use the
	\eqn{(2s+1)^{-2}} (this was the default in versions prior to
	1.7.4)
}
  \item{\dots}{
	Other arguments
}
}
\details{1. If \code{binwidth=0} the function first computes the
	`best' linear running mean binwidth (bandwidth)
	for the smooth of the localized autocovariance.
	2. The function computes the localized autocovariance
	smoothed with a running mean with the selected binwidth.
	Then, the variance of \eqn{\hat{c}(z, \tau)} is
	computed for the selected value of time z=nz and for the
	lags specified (in \code{var.lag.max}). The results are
	returned in an object of class \code{lacfCI}.

	Note, this function
	computes and plots localized autocovariances for a particular
	and fixed time location. Various other plots, including
	perspective plots or the localized autocovariance function
	over all time can be found in the \code{costat} package.
	(Indeed, this function returns a \code{lacfCI} object that
	contains a \code{lacf} object, and interesting plots
	can be plotted using the \code{plot.lacf} function within
	\code{costast}.
}
\value{An object of class \code{lacfCI}. This is a list with
	the following components.
\item{lag}{The lags for which the localized autocovariance variance is
	computed}
\item{cvar}{The variances associated with each localized autocovariance}
\item{the.lacf}{The \code{lacf} class object that contains the
	localized autocovariances themselves. This object
	can be handled/plotted/etc using the functions
	in the \code{costat} package although \code{plot.lacfCI}
	contains much of the functionality of \code{plot.lacf}.} 
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

\seealso{\code{\link{plot.lacfCI}}, \code{\link{print.lacfCI}},
\code{\link{summary.lacfCI}}
}
\examples{
#
# Do localized autocovariance on a iid Gaussian sequence
#
tmp <- Rvarlacf(rnorm(256), nz=125)
#
# Plot the localized autocovariances over time (default plot, doesn't
# produce CIs)
#
\dontrun{plot(tmp)}
#
# You should get a plot where the lag 0 acs are all near 1 and all the
# others are near zero, the acfs over time. 
#
\dontrun{plot(tmp, plotcor=FALSE, type="acf")}
#
# This plots the autocovariances (note: \code{plotcor=FALSE}) and the
# type of plot is \code{"acf"} which means like a regular ACF plot, except
# this is c(125, tau), ie the acf localized to time=125 The confidence
# intervals are also plotted. 
# The plot subtitle indicates that it is c(125, tau) that is being plotted
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
