\name{plot.lacfCI}
\alias{plot.lacfCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot confidence intervals for localized autocovariance for
	locally stationary time series.
}
\description{Plot the localized autocovariance and
	approximate confidence intervals.
}
\usage{
\method{plot}{lacfCI}(x, plotcor = TRUE, type = "line",
    lags = 0:as.integer(10 * log10(nrow(x$lacf))), tcex = 1,
    lcol = 1, llty = 1, ylim = NULL, segwid = 1,
    segandcross = TRUE, conf.level = 0.95, plot.it = TRUE,
    xlab, ylab, sub, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The \code{lacfCI} object you wish to plot, e.g. produced
	by the \code{Rvarlacf} function.
} 
  \item{plotcor}{If \code{TRUE} then autocorrelations are plotted,
	if \code{FALSE} then autocovariances are. Note: not all combinations
	of types of plot and \code{plotcor} are valid, but many are.
}
  \item{type}{This can be one of three values \code{"line"},
	\code{"persp"} or \code{"acf"}. The value
	\code{"acf"} produces a plot like the regular \code{acf}
	function, but note, the values plotted are from a localized
	autocovariance function centred at the time location
	contained in the \code{object} object (and that time appears in
	the subtitle). This is the only plot that also plots the
	confidence intervals. The \code{"line"} plot plots 
	autocorrelations (only) for the specified \code{lags}
	and does this over all time for the whole extent of the
	series. This plot is useful to see if the autocorrelations
	are changing over time. The final option,
	\code{"persp"} produces a perspective plot of the autocovariance
	or autocorrelations. Arguments can be supplied (\code{theta},
	\code{phi}) to rotate the perspective plot, as it can be
	sometimes hard to visualize the plot.
}
  \item{lags}{The lags that you wish to display. This
	should be a list of non-negative integers, but not necessarily
	consecutive.
}
  \item{tcex}{On the \code{"line"} plot this argument controls
	the expansion of the font for the labels on the lines.
	So, setting \code{tcex=2}, for example, will double the size
	of these. These labels visually indicate which line corresponds
	to which lag.
}
  \item{lcol}{On the \code{"line"} plot, this argument
	controls the colour of the lines that are used to show
	the acfs.
}
  \item{llty}{As \code{lcol} but for line types.
}
  \item{ylim}{The vertical limits of the plot.
}
  \item{segwid}{On the \code{"acf"} plot, this argument controls
	the widths of the little acf segments that connect the x-axis
	with the acf values.
}
  \item{segandcross}{If \code{TRUE} then a small diamond is plotted
	at the location of the acf, to make it clearer.
}
  \item{conf.level}{The confidence level of the confidence intervals.
}
  \item{plot.it}{If \code{FALSE} then no plot is produced. This
	can be used if you merely want to extract the relevant
	acf values (which are returned).
}
  \item{xlab}{X-axis label, constructed internally if not supplied}
  \item{ylab}{Y-axis label, constructed internally if not supplied}
  \item{sub}{A subtitle for the plot}
  \item{\dots}{Other arguments to the main \code{plot} command.
}
}
\details{This function can plot the localized autocovariance in
	three ways. Like a regular acf plot (but obviously a slice
	out of a time-varying autocovariance, not the regular acf),
	a line plot which shows the acfs over time and a perspective
	plot which can plot the estimate of \eqn{c(z, \tau)} as a
	2D function. Currently, the confidence intervals can only
	be displayed on the \code{"acf"} type plot.
}
\value{A vector of the extracted acfvals invisibly returned.
	Note: what is returned depends on the arguments, what is
	returned is what would have been plotted if \code{plot.it}
	were \code{TRUE}
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

\seealso{
\code{\link{Rvarlacf}}
}
\examples{
#
# Simulate a TVAR(1) process
#
x <- tvar1sim()
#
# Computes its time-localized autocovariance and confidence intervals
# Note: smoothing is done automatically!
#
x.lacf <- Rvarlacf(x=x, nz=50, var.lag.max=20)
#
# Now plot this, plot covariances as an acf plot, with the CIs
#
\dontrun{plot(x.lacf, type="acf", plotcor=FALSE)}
#
# Now plot it as a line plot, as correlations and can't do CIs
#
\dontrun{plot(x.lacf)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
