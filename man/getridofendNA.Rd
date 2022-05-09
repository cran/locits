\name{getridofendNA}
\alias{getridofendNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Replaces all NAs in vector by 0
}
\description{Replaces all NAs in vector by 0
}
\usage{
getridofendNA(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	Vector that might contain NAs
}
}
\details{Originally, this function did something more complex,
but now it merely replaces NAs by 0
}
\value{The same vector as \code{x} but with NAs replaced by 0
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

\seealso{\code{\link{HwdS}}
}
\examples{
#
#
#
x <- c(3, 4, 6, NA, 3)
getridofendNA(x)
#[1] 3 4 6 0 3
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
