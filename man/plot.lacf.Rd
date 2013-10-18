\name{plot.lacf}
\alias{plot.lacf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot localized autocovariance (lacf) object.
}
\description{Produces various ways of looking at a localized autocovariance
	(lacf) object.
}
\usage{
\method{plot}{lacf}(x, plotcor = TRUE, type = "line",
        lags = 0:min(as.integer(10 * log10(nrow(x$lacf))), ncol(x$lacf) - 1),
        tcex = 1, lcol = 1, llty = 1, the.time = NULL, plot.it=TRUE,
	xlab, ylab, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The localized autocovariance object you want to plot (lacf)
}
  \item{plotcor}{If TRUE then plot autocorrelations, otherwise plot
	autocovariances.
}
  \item{type}{The lacf objects are fairly complex and so there are
	different ways you can plot them. The \code{type}s are
	\code{line}, \code{persp} or \code{acf}, see the details for
	description. Note that the \code{line} plot only works with
	correlations currently.
}
  \item{lags}{The \code{lags} that you wish included in the plot.
	The default is all the lags from 0 up to the maximum that
	is used in the R acf plot 
}
  \item{tcex}{In the \code{line} plot lines are plotted that indicate
	the time-varying correlation. Each lag gets a different line
	and the lines are differentiated by the lag id being placed
	at intervals along the line. This argument changes the
	size of those ids (numbers).
}
  \item{lcol}{Controls the colours of the lines in the \code{line}
	plot.
}
  \item{llty}{Controls the line types of the lines in the \code{line}
	plot.
}
  \item{the.time}{If the \code{acf} plot is chosen then you have to specify
	a time point about which to plot the acf. I.e. in general this
	funcion's lacf argument is a 2D function: \eqn{c(t, \tau)}, the
	\code{acf} plot produces a plot like the regular acf function
	and so you have to turn the 2D \eqn{c(t, \tau)} into a
	1D function \eqn{c(t_0, \tau)} by specifying a fixed time
	point \eqn{t_0}.
}
  \item{plot.it}{If \code{TRUE} the plot is produced and displayed. If
	\code{FALSE} then no plot is produced but the autocovariance
	or autocorrelation values that would have been produced are
	returned as numerical values instead. This means that
	this function is an extractor function for the \code{lacf}
	class object.
}
  \item{xlab}{X-axis label, constructed internally if not supplied}
  \item{ylab}{Y-axis label, constructed internally if not supplied}
  \item{\dots}{Other arguments to plot.
}
}
\details{This function produces pictures of the
	two-dimensional time-varying autocovariance
	or autocorrelation, \eqn{c(t, \tau)},
	of a locally stationary time series.
	There are three types of plot depending on the argument to
	the \code{type} argument.

	The \code{line} plot draws the autocorrelations as a series of
	lines, one for each lag, as lines over time. E.g. a sequence
	#of lines \eqn{c(t, \tau_i)} is drawn, one for each \eqn{\tau_i}.
	The zeroth lag line is the autocorrelation at lag 0 which is
	always 1. By default all the lags are drawn which can result
	in a confusing picture. Often, one is only interested in the low
	level lags, so only these can be plotted by changing the \code{lags}
	argument and any selection of lags can be plotted. The colour
	and line type of the plotted lines can be changed with the
	\code{lcol} and the \code{llty} arguments.

	The \code{acf} plot produces pictures similar to the standard
	R \code{acf()} function plot. However, the regular acf is a
	1D function, since it is defined to be constant over all time.
	The time-varying acf supplied to this function is not constant
	over all time (except for stationary processes, theoretically).
	So, this type of plot requires the user to specify a fixed
	time at which to produce the plot, and this is supplied by
	the \code{the.time} argument.

	The \code{persp} plot plots the 2D function \eqn{c(t, \tau)}
	as a perspective plot. 
}
\value{For the \code{acf} type plot the acf values are returned
	invisibly. For the other types nothing is returned.
}
\references{
Cardinali, A. and Nason, G.P. (2012) Costationarity of Locally 
Stationary Time Series using \code{costat}.

Cardinali, A. and Nason, G.P. (2010) Costationarity of locally stationary
time series. \emph{J. Time Series Econometrics}, \bold{2}, Issue 2, Article 1.

Nason, G.P. (2013) A test for second-order stationarity and
approximate confidence intervals for localized autocovariances
for locally stationary time series. \emph{J. R. Statist. Soc.} B,
\bold{75}, 879-904.
}
\author{
G.P. Nason
}

\seealso{\code{\link{lacf}}}
\examples{
#
# Make some dummy data, e.g. white noise
#
v <- rnorm(256)
#
# Compute the localized autocovariance (ok, the input is stationary
# but this is just an example. More interesting things could be achieved
# by putting the results of simulating from a LSW process, or piecewise
# stationary by concatenating different stationary realizations, etc.
#
vlacf <- lacf(v, lag.max=30)
#
# Now let's do some plotting of the localized autocovariance
#
\dontrun{plot(vlacf, lags=0:6)}
#
# Should get a plot where lag 0 is all up at value 1, and all other
# autocorrelations are near zero (since its white noise).
#
#
# How about just looking at lags 0, 2 and 4, and some different colours.
#
\dontrun{plot(vlacf, lags=c(0,2,4), lcol=c(1,2,3))}
#
# O.k. Let's concentrate on time t=200, let's look at a standard acf
# plot near there.
#
\dontrun{plot(vlacf, type="acf", the.time=200)}
#
# Now plot the autocovariance, rather than the autocorrelation.
#
\dontrun{plot(vlacf, type="acf", the.time=200, plotcor=FALSE)}
#
# Actually, the plot doesn't look a lot different as the series is white
# noise, but it is different if you look closely.
}
\keyword{ts}
