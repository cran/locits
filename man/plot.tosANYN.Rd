\name{plot.tosANYN}
\alias{plot.tosANYN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produces a graphical representation of the results of a
	test of stationarity from a \code{tosANYN} object.
}
\description{After a test of stationarity (e.g. \code{\link{hwtos}})
	is applied to a time series it generates a results object of
	class \code{tosANYN}. This function takes objects of that
	class and produces a graphical representation of the test.
}
\usage{
\method{plot}{tosANYN}(x, sub = NULL, xlab = "Time",
    arrow.length = 0.05, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The \code{tosANYN} class object, the results of the
	test of stationarity that you wish to plot.
}
  \item{sub}{An argument to change the subtitle.
}
  \item{xlab}{An argument to change the x-axis label.
}
  \item{arrow.length}{The length of the edges of the arrow
	head (in inches). Note that this is the argument that is
	supplied as the \code{length} argument of the \code{arrow}
	function that is called by this routine to draw the arrows.
}
  \item{verbose}{If \code{TRUE} then some meaningless debugging
	information is printed.
}
  \item{\dots}{
	Other arguments to the main \code{ts.plot} routine that
	does the plotting.
}
}
\details{The following things are usually plotted. 1. The time series
	that was investigated. The left-hand axes is that for the time series.
	The horizontal axis is time (but just integers indexing).
	If the series was deemed stationary by the test then that's it
	except that the subtitle indicates that no Haar wavelet coefficients
	were rejected as being nonzero.

	If the test indicated that the series was nonstationary then
	the subtitle indicates this by stating the number of rejections.
	Then graphical
	representations of any significant Haar wavelet coefficients
	are plotted as double-headed red horizontal arrows on the plot.
	The horizontal extent corresponds to the support of the underlying
	wavelet. The vertical position of the arrows gives an indication
	of the wavelet periodogram scale where the significant coefficient
	was found. The wavelet periodogram scales are indexed by the right
	hand axis, and beware, the numbers might not be consecutive, but
	the will be ordered (so e.g. if no signficant coefficients were
	discovered at wavelet periodogram scale level 6, then that scale/axis
	label will not appear). The scale within the Haar wavelet transform
	is indicated by the vertical position WITHIN ticks between
	wavelet periodogram scales (ie, there are TWO scales: the wavelet
	periodogram scale that is currently being analyzed, and the
	Haar wavelet transform scale within the periodogram scale).
	So, if two right hand axis labels are, e.g., 4 and 5, and
	horizontal arrows appear between these two they actually correspond
	to different Haar wavelet transform scales AT wavelet periodogram
	level 4. It is not usually possible to tell precisely which
	Haar wavelet transform scale the coefficients can come from,
	but the information can be extracted from the \code{\link{summary.tosANYN}}
	function which lists this.
}
\value{None.
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

\seealso{\code{\link{hwtos}}, \code{\link{summary.tosANYN}}
}
\examples{
#
# Produces an interesting plot with high probability
#
#
# Note that the input time series is two concatenated white noise
# sequences with very different variances.
#
\dontrun{answer <- hwtos(c(rnorm(256), rnorm(256, sd=5)))}
\dontrun{plot(answer)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
