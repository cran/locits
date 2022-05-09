\name{idlastzero}
\alias{idlastzero}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Return the index of the last zero in a vector
}
\description{Return the index of the last zero in a vector, otherwise
	stop and return errror message. A helper routine for
	\code{\link{mkcoef}}.
}
\usage{
idlastzero(v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{Vector you wish to investigate
}
}
\value{The index within \code{v} of the last (right-most or one with
	the largest index) zero.
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
\code{\link{mkcoef}}
}
\examples{
idlastzero(c(3,4,5,0,9))
#[1] 4
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
