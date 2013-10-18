\name{ewspec3}
\alias{ewspec3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute evolutionary wavelet spectrum of a time series.
}
\description{This function is a development of the \code{ewspec}
	function from \code{wavethresh} but with more features.
	The two new features are: the addition of running mean smoothing
	and autoreflection which mitigates the problems caused in
	\code{ewspec} which performed periodic transforms on
	data (time series) which were generally not periodic. 
}
\usage{
ewspec3(x, filter.number = 10, family = "DaubLeAsymm",
    UseLocalSpec = TRUE, DoSWT = TRUE, WPsmooth = TRUE,
    WPsmooth.type = "RM", binwidth = 5, verbose = FALSE,
    smooth.filter.number = 10, smooth.family = "DaubLeAsymm",
    smooth.levels = 3:WPwst$nlevels - 1, smooth.dev = madmad,
    smooth.policy = "LSuniversal", smooth.value = 0,
    smooth.by.level = FALSE, smooth.type = "soft",
    smooth.verbose = FALSE, smooth.cvtol = 0.01,
    smooth.cvnorm = l2norm, smooth.transform = I,
    smooth.inverse = I, AutoReflect = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The time series you want to compute the evolutionary
	wavelet spectrum for.
}
  \item{filter.number}{Wavelet filter number underlying the analysis
	of the spectrum (see \code{filter.select} or \code{wd} for more
	details).
}
  \item{family}{Wavelet family. Again, see \code{filter.select} or \code{wd}
	for more details.
}
  \item{UseLocalSpec}{ As \code{ewspec}, should usually leave as is.
}
  \item{DoSWT}{As \code{ewspec}, should usually leave as is.
}
  \item{WPsmooth}{If \code{TRUE} then smoothing is applied to
	the wavelet periodogram (and hence spectrum).
}
  \item{WPsmooth.type}{The type of periodogram smoothing. 
	If this argument is \code{"RM"} then running mean
	linear smoothing is used.
	Otherwise, wavelet shrinkage as in \code{ewspec} is
	used.
}
  \item{binwidth}{If the periodogram smoothing is \code{"RM"} then
	the this argument supplies the \code{binwidth} or number
	of consecutive observations used in the running mean smooth.
}
  \item{verbose}{If \code{TRUE} then messages are produced. If 
	\code{FALSE} then they are not.
}
  \item{smooth.filter.number}{If wavelet smoothing of the wavelet
	periodogram is used then this specifies the index number of
	wavelet to use, exactly as \code{ewspec}.
}
  \item{smooth.family}{If wavelet smoothing of the wavelet
        periodogram is used then this specifies the family of
        wavelet to use, exactly as \code{ewspec}.
}
  \item{smooth.levels}{If wavelet smoothing of the wavelet
        periodogram is used then this specifies the levels to
	smooth, exactly as \code{ewspec}.
}
  \item{smooth.dev}{If wavelet smoothing of the wavelet
        periodogram is used then this specifies  deviance used
	to compute smoothing thresholds, exactly as \code{ewspec}.
}
  \item{smooth.policy}{If wavelet smoothing of the wavelet
        periodogram is used then this specifies the policy
        of wavelet shrinkage to use, exactly as \code{ewspec}.
}
  \item{smooth.value}{If wavelet smoothing of the wavelet
        periodogram is used then this specifies the value of the
	smoothing parameter for some policies, exactly as \code{ewspec}.
}
  \item{smooth.by.level}{If wavelet smoothing of the wavelet
        periodogram is used then this specifies whether level-by-level
	thresholding is applied, or one threshold is applied to
	all levels, exactly as \code{ewspec}.
}
  \item{smooth.type}{If wavelet smoothing of the wavelet
        periodogram is used then this specifies the type of
	thresholding, "hard" or "soft", exactly as \code{ewspec}.
}
  \item{smooth.verbose}{If wavelet smoothing of the wavelet
        periodogram is used then this specifies whether or not
	verbose messages are produced during the smoothing,
        exactly as \code{ewspec}.
}
  \item{smooth.cvtol}{If wavelet smoothing of the wavelet
        periodogram is used then this specifies a tolerance
	for the cross-validation algorithm if it is specified
	in the \code{smooth.policy}, exactly as \code{ewspec}.
}
  \item{smooth.cvnorm}{Ditto to the previous argument, but this
	one supplies the norm used by the cross-validation.
}
  \item{smooth.transform}{If wavelet smoothing of the wavelet
        periodogram is used then this specifies whether a transform
	is used to transform the periodogram before smoothing, 
        exactly as \code{ewspec}.
}
  \item{smooth.inverse}{Should be the mathematical inverse of
	the \code{smooth.transform} argument.
}
  \item{AutoReflect}{Whether the series is internally reflected before
	application of the wavelet transforms. So, \code{x} becomes
	\code{c(x, rev(x))} which is a periodic sequence. After
	estimation of the spectrum the second-half of the spectral
	estimate is junked (because it is a reflection of the first
	half). However, the estimate is better. This argument improves
	over \code{ewspec} where poor estimates near boundaries were
	obtained because the transforms assume periodicity but most
	time series are not (and X_1 and X_T are very different, etc).
}
}
\value{Precisely the same kind of output as \code{ewspec}.
}


\references{Nason, G.P. (2013) A test for second-order stationarity and
	approximate confidence intervals for localized autocovariances
	for locally stationary time series. \emph{J. R. Statist. Soc.} B,
	\bold{75}, 879-904.}
\author{
Guy Nason.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{AutoBestBW}}, \code{\link{lacf}}
}
\examples{
#
# Generate time series
#
x <- tvar1sim()
#
# Compute its evolutionary wavelet spectrum, with linear running mean smooth
#
x.ewspec3 <- ewspec3(x)
#
# Plot the answer, probably its a bit variable, because the default bandwidth
# is 5, which is probably inappropriate for many series
#
\dontrun{plot(x.ewspec3$S)}
#
# Try a larger bandwidth
#
x.ewspec3 <- ewspec3(x, binwidth=100)
#
# Plot the answer, should look a lot smoother
#
# Note, a lot of high frequency power on the right hand side of the plot,
# which is expected as process looks like AR(1) with param of -0.9
#
\dontrun{plot(x.ewspec3$S)}
#
# Do smoothing like ewspec (but additionally AutoReflect)
#
x.ewspec3 <- ewspec3(x, WPsmooth.type="wavelet")
#
# Plot the results
#
\dontrun{plot(x.ewspec3$S)}
#
# Another possibility is to use AutoBestBW which tries to find the best
# linear smooth closest to a wavelet smooth. This makes use of ewspec3
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
