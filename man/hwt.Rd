\name{hwt}
\alias{hwt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute a Haar wavelet transform for data of arbitrary n length
}
\description{Function computes Haar wavelet and scaling function coefficients
	for data set of any length. Algorithm computes every possible
	coefficient that it can for both decimated and nondecimated
	versions of the transform.
}
\usage{
hwt(x, type = c("wavelet", "station"), reindex = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of length n, where n is a positive integer.
	This is the data that you wish to compute the Haar wavelet
	transform for.
}
  \item{type}{The type of transform, either the decimated or nondecimated
	algorithm.
}
  \item{reindex}{If TRUE then the routine attempts to match scales
	with the usual dyadic transform, \code{wd}. If FALSE then the
	coefficients that are returned are "as is"
}
}
\details{Essentially, this algorithm attempts to compute every
	possible Haar wavelet coefficient. For example, if the length
	of the input series was 6 then this means that three coefficients
	at the finest scale can be computed using the first, second and
	third pair of input data points using the weights
	\code{c(1, -1)/sqrt(2)}. However, from the three coefficients
	that result from this, there is only one pair, so only one
	"next coarser" coefficient can be computed.

	The \code{reindex} option is subtle. Essentially, it tries to
	ensure that the returned coefficients end up at the same
	scales as if a data set of the next highest dyadic length was
	analyzed by the wd function. E.g. if the length of the
	series was 10 then with \code{reindex=FALSE} (default)
	only three levels are returned for each of the wavelet and
	scaling coefficients. If \code{reindex=TRUE} then the number
	of levels returned would be as if \code{wd} analysed a data set
	of length 16 (the smallest dyadic number larger than 10).
	The wd levels would be zero to three and this is what
	would be returned in this function if \code{reindex=TRUE}.
	However, note, in this case, the coarsest level coefficient
	happens to be NULL (or not computable). One can view the algorithm
	as computing a partial transform of 10 of the 16 elements and
	substituting NA for anything it can't compute. 
}
\value{An object of class hwtANYN which is a list with the following
	components.
\item{c}{The scaling function coefficients. This is a list of length
	\code{nlevels} which contains the scaling function coefficients.
	The coarsest scale coefficients are to be found in the lowest-indexed
	slots of the list (e.g. c[[1]]) and increasing slot index corresponds
	to finer scales. So, \code{c[[length(c)]]} corresponds to the finest
	coefficients. Note, an entry in the slot can also be NULL. This
	indicates that no coefficients could be calculated at this scale,
	usually the coarsest.}
\item{d}{A for \code{c} but for wavelet coefficients.}
\item{nlevels}{The number of scale levels in the Haar wavelet decomposition.
	if \code{reindex=TRUE} then this number will be the log to base 2
	of the smallest power of two larger than the length of the input
	vector \code{x}.}
\item{type}{Whether a decimated wavelet transform has been computed
	(\code{"wavelet"}) or a nondecimated transform (\code{"station"}).
	Note, the name of the argument \code{"station"} has been chosen
	to coincide with the \code{type} in the regular wavelet
	transform computed by \code{wd}.}
\item{reindex}{Either \code{TRUE} or \code{FALSE}. If \code{TRUE} then
	the scale levels correspond directly to those computed by
	\code{wd}, the regular wavelet transform. If \code{FALSE} then
	the levels returned in \code{c} and \code{d} are just indexed from
	the first non-null level.}
}

\references{
Nason, G.P. (2013) A test for second-order stationarity and
        approximate confidence intervals for localized autocovariances
        for locally stationary time series. \emph{J. R. Statist. Soc.} B,
        \bold{75}, 879-904.
	\doi{10.1111/rssb.12015}

Priestley, M.B. and Subba Rao (1969) A test for non-stationarity of
        time series. \emph{J. R. Statist. Soc. B}, \bold{31}, 140-149.

von Sachs, R. and Neumann, M.H. (2000) A wavelet-based test for
        stationarity. \emph{J. Time Ser. Anal.}, \bold{21},
        597-613.
}

\author{
G. P. Nason
}

\seealso{ \code{\link{hwtos}},
	\code{\link{plot.hwtANYN}}, \code{\link{print.hwtANYN}},
	\code{\link{summary.hwtANYN}}
}
\examples{
#
# Generate test data set of length 5 (note, NOT a power of two)
#
v2 <- rnorm(5)
#
# Compute its Haar transform
#
v2hwt <- hwt(v2)
#
# How many levels does it have?
#
nlevelsWT(v2hwt)
#
# What are the coarsest scale wavelet coefficients
#
v2hwt$d[[1]]
#
# What are the finest scale scaling function coefficients
#
v2hwt$c[[nlevels(v2hwt$c)-1]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
\keyword{smooth}
\keyword{nonlinear}
