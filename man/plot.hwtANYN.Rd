\name{plot.hwtANYN}
\alias{plot.hwtANYN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots the transform contained in an \code{hwtANYN} object.
}
\description{An \code{hwtANYN} object contains the results of a Haar
	wavelet transform computed on an object of non-dyadic length.
	It is the equivalent of the \code{wd} object for non-dyadic
	vectors for Haar wavelets. Note, the plot can only be carried
	out where the \code{reindex} slot of the object is \code{TRUE}.
}
\usage{
\method{plot}{hwtANYN}(x, xlabvals, xlabchars, ylabchars, first.level = 1,
    main = "Haar Wavelet Coefficients", scaling = c("global", "by.level"),
    rhlab = FALSE, sub, NotPlotVal = 0.005, xlab = "Translate",
    ylab = "wd-equivalent Resolution Level", miss.coef.col = 2,
    miss.coef.cex = 0.5, miss.coef.pch = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The \code{hwtANYN} object containing the Haar wavelet transform
	coefficients you wish to plot.
}
  \item{xlabvals}{Coordinates of x-axis labels you wish to add.
}
  \item{xlabchars}{Labels to be printed at the x-axis labels specified.
}
  \item{ylabchars}{Y-axis labels}
  \item{first.level}{Specifies the coarsest level to be plotted.
}
  \item{main}{Specify a different main title for the plot.
}
  \item{scaling}{How coefficients will be scaled on the plot.
	This can be two arguments \code{"global"} where all coefficients
	are plotted to the same scale and \code{"by.level"} where all
	coefficients on the same resolution level are plotted to the same
	scale, but coefficients on different resolution levels might be
	of different scales. 
}
  \item{rhlab}{If \code{TRUE} then the scale factor used for each level is
	shown.
}
  \item{sub}{Specify a different subtitle for the plot.
}
  \item{NotPlotVal}{Coefficients will not be plotted if their scaled
	height is less than \code{NotPlotVal} in absolute value. This is a
	useful way to completely suppress very small coefficient values.
}
  \item{xlab}{Specify the x-axis label.
}
  \item{ylab}{Specify the y-axis label.
}
  \item{miss.coef.col}{What color to plot "missing coefficients" in.
}
  \item{miss.coef.cex}{How big to plot the "missing coefficients" symbol.
}
  \item{miss.coef.pch}{The type of plotting character used to plot the
	"missing coefficients".
}
  \item{\dots}{Other arguments to plot.
}
}
\details{A plot of the different wavelet coefficients at the scales
	ranging from \code{first.level} to the finest scale. Note, in this
	plot the coefficients are NOT aligned with time at different
	scales in the same way as in the \code{wd} type plot
	- except the finest scale.

	The Haar wavelet transform objects that this function plots
	are obtained originally from vectors of non-dyadic length. One
	can think of such a vector as a sub-vector of a longer vector
	of dyadic length. E.g. if your vector is of length 35 then it is
	a sub-vector of a vector of 64 (the next highest power of two).
	So, you can think of the Haar wavelet transform being of a vector
	of length 64 where 64-35=29 of the observations are missing. These
	missing observations "contribute" to wavelet (and scaling function)
	coefficients that are missing. This function has the ability to
	plot the "missing" coefficients, by default as small red triangles.
	The user can control the colour, size and plotting character of the
	missing observations.
}
\value{A single vector of length the number of levels plotted containing
	the value of the maximum absolute coefficient value.
}
\references{
Nason, G.P. (2013) A test for second-order stationarity and
        approximate confidence intervals for localized autocovariances
        for locally stationary time series. \emph{J. R. Statist. Soc.} B,
        \bold{75}, 879-904.

Priestley, M.B. and Subba Rao (1969) A test for non-stationarity of
        time series. \emph{J. R. Statist. Soc. B}, \bold{31}, 140-149.

von Sachs, R. and Neumann, M.H. (2000) A wavelet-based test for
        stationarity. \emph{J. Time Ser. Anal.}, \bold{21},
        597-613.
}
\author{
G. P. Nason
}


\seealso{\code{\link{hwt}}, \code{\link{print.hwtANYN}}
}
\examples{
#
# Generate test data of length 82
#
v3 <- rnorm(82)
#
# Compute Haar wavelet transform, note reindex has to be true for subsequent
# plot.
#
v3.hwt <- hwt(v3, reindex=TRUE)
#
#
\dontrun{plot(v3.hwt)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{smooth}
