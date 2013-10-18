\name{Cvarip2}
\alias{Cvarip2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes variance of Haar wavelet coefficients of
	wavelet periodogram using C code.
}
\description{Performs precisely the same role as \code{\link{varip2}}
	except it is implemented internally using C code and hence
	is much faster.
}
\usage{
Cvarip2(i, p, ll, S, Pmat, PsiJL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{Scale parameter of Haar wavelet analyzing periodogram.
	Scale 1 is the finest scale.
}
  \item{p}{Location parameter of Haar wavelet analyzing periodogram
}
  \item{ll}{Scale of the raw wavelet periodogram being analyzed.
}
  \item{S}{Estimate of the spectrum, under the assumption of stationarity.
        So, this is just a vector of (possibly) J scales (which is often
        the usual spectral estimate averaged over time). Note: that the
        main calling function, \code{\link{hwtos2}}, actually passes
        \code{maxD} levels.
}
  \item{Pmat}{Matrix version of autocorrelation wavelet computed
	using the \code{PsiJmat} function in \code{wavethresh}
}
  \item{PsiJL}{True length of the autocorrelation wavelets
	in the \code{Pmat} matrix. This can be obtained simply
	by using the list version of the ac wavelet (computed
	by \code{PsiJ}) and applying \code{sapply}.
}
}
\value{The list returned from the \code{.C} calling function.
The only object of real interest is the \code{ans} component
which contains the variance.
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
\code{\link{varip2}}
}
\examples{
#
# See example from varip2
#
#
my.Pmat <- PsiJmat(-5, filter.number=1, family="DaubExPhase")
my.PsiJ <- PsiJ(-5, filter.number=1, family="DaubExPhase")
my.PsiJL <- sapply(my.PsiJ, "length")
Cvarip2(i=1, p=10, ll=2, S=c(1/2,1/4,1/8,1/16,1/32),
    Pmat=my.Pmat, PsiJL=my.PsiJL)
#
# Gives answer 1.865244, which is the same as given in the example for varip2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
