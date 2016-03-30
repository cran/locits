\name{EstBetaCov}
\alias{EstBetaCov}
\title{Compute estimate of wavelet periodogram and the estimate's
        covariance matrix.
}
\description{An estimate of the wavelet periodogram at a location
        \code{nz} is generated. This is obtained by first computing the
        empirical raw wavelet periodogram by squaring the results of
        the nondecimated wavelet transform of the time series. Then
        a running mean smooth is applied.  }
\usage{
EstBetaCov(x, nz, filter.number = 1, family = "DaubExPhase", smooth.dev = var,
    AutoReflect = TRUE, WPsmooth.type = "RM", binwidth = 0, mkcoefOBJ,
    ThePsiJ, Cverbose = 0, verbose = 0, OPLENGTH = 10^5, ABB.tol = 0.1,
    ABB.plot.it = FALSE, ABB.verbose = 0, ABB.maxits = 10, do.init = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The time series for which you wish to have the estimate for.
}
  \item{nz}{The time point at which you want the estimate computed at.
	This is an integer ranging from one up to the length of the time
	series.
}
  \item{filter.number}{The analysis wavelet (the wavelet periodogram is
	computed using this to form the nondecimated wavelet coefficients)
}
  \item{family}{The family of the analysis wavelet.
}
  \item{smooth.dev}{The deviance function used in smoothing via the
	internal call to the \code{\link{ewspec3}} function.
}
  \item{AutoReflect}{Whether better smoothing is to be obtained by
	AutoReflection to mitigate the effects of using periodic
	transforms on non-periodic data. See \code{\link{ewspec3}}
}
  \item{WPsmooth.type}{The type of wavelet periodogram smoothing. For here
	leave the option at \code{"RM"} otherwise unpredictable results can occur
}
  \item{binwidth}{The running mean length. If zero then a good bandwidth
	will be chosen using the \code{\link{AutoBestBW}} function.
}
  \item{mkcoefOBJ}{If this argument is missing then it is computed internally
	using the \code{\link{mkcoef}} function which computes discrete
	wavelets. If this function is going to be repeatedly called then
	it is more efficient to supply this function with a precomputed
	version.
}
  \item{ThePsiJ}{As for \code{mkcoefOBJ} argument but for the
	autocorrelation wavelet and the function \code{\link{PsiJ}}.
}
  \item{Cverbose}{This function called the C routine \code{CstarIcov}
	if you set \code{Cverbose} to true then the routine instructs
	the C code to produce debugging messages.
}
  \item{verbose}{If TRUE then debugging messages from the R code are
	produced.
}
  \item{OPLENGTH}{Subsidiary parameters for potential call to
	\code{\link{PsiJ}} function
}
  \item{ABB.tol}{Tolerance to be passed to \code{\link{AutoBestBW}} function.
}
  \item{ABB.plot.it}{Argument to be passed to \code{\link{AutoBestBW}}
	\code{plot.it} argument.
}
  \item{ABB.verbose}{Argument to be passed to \code{\link{AutoBestBW}}
	\code{verbose} argument.
}
  \item{ABB.maxits}{Argument to be passed to \code{\link{AutoBestBW}}
	\code{maxits} argument.
}
  \item{do.init}{Initialize stored statistics, for cache hit rate info.
}
  \item{\dots}{Other arguments that are passed to the
	\code{\link{ewspec3}} function.
}
}
\details{First optionally computes a good bandwidth using the
	\code{\link{AutoBestBW}} function. Then
	computes raw wavelet periodogram using \code{\link{ewspec3}}
	using running mean smoothing with the \code{binwidth} bandwith
	(which might be automatically chosen). This computes the estimate
	of the wavelet periodogram at time \code{nz}. The covariance matrix
	of this estimate is then computed in C using the
	\code{CstarIcov} function and this is returned.
}
\value{A list with two components:
	\item{betahat}{A vector of length \eqn{J} (the number of scales
	in the wavelet periodogram, which is \eqn{\log_2}
	of the number of observations \eqn{T}}
	\item{Sigma}{A matrix of dimensions \eqn{J\times J}
	which is the covariance
	of \eqn{\hat{\beta}_j} with \eqn{\hat{\beta}_\ell}.}
}
\references{
Nason, G.P. (2013) A test for second-order stationarity and
        approximate confidence intervals for localized autocovariances
        for locally stationary time series. \emph{J. R. Statist. Soc.} B,
        \bold{75}, 879-904.
}
\author{Guy Nason
}

\seealso{
\code{\link{AutoBestBW}}, \code{\link{ewspec3}}
}
\examples{
#
# Small example, not too computationally demanding on white noise
#
myb <- EstBetaCov(rnorm(64), nz=32)
#
# Let's see the results (of my run)
#
\dontrun{myb$betahat}
#[1] 0.8323344 0.7926963 0.7272328 1.3459313 2.1873395 0.8364632
#
# For white noise, these values should be 1 (they're estimates)
\dontrun{myb$Sigma}
#            [,1]        [,2]        [,3]       [,4]        [,5]        [,6]
#[1,] 0.039355673 0.022886994 0.008980497 0.01146325 0.003211176 0.001064377
#[2,] 0.022886994 0.054363333 0.035228164 0.06519112 0.017146883 0.006079162
#[3,] 0.008980497 0.035228164 0.161340373 0.38326812 0.111068916 0.040068318
#[4,] 0.011463247 0.065191118 0.383268115 1.31229598 0.632725858 0.228574601
#[5,] 0.003211176 0.017146883 0.111068916 0.63272586 1.587765187 0.919247252
#[6,] 0.001064377 0.006079162 0.040068318 0.22857460 0.919247252 2.767615374
#
# Here's an example for T (length of series) bigger, T=1024
#
\dontrun{myb <- EstBetaCov(rnorm(1024), nz=512)}
#
# Let's look at results
#
\dontrun{myb$betahat}
# [1] 1.0276157 1.0626069 0.9138419 1.1275545 1.4161028 0.9147333 1.1935089
# [8] 0.6598547 1.1355896 2.3374615
#
# These values (especially for finer scales) are closer to 1
#
}
\keyword{ts}
