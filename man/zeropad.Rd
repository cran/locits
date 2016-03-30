\name{zeropad}
\alias{zeropad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Intersperse zeroes in a vector.
}
\description{Take a vector of any length and then intersperse zeroes
	between existing elements (and add a further zero to the end).
}
\usage{
zeropad(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector that you want to intersperse zeros into.
}
}
\details{Title says it all.
}
\value{A vector whose odd elements are just \code{x} and whose
	even elements are zeroes.
}
\author{
G.P. Nason
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{hwt}}
}
\examples{
#
# Operate on a test set
#
v <- 1:3
zeropad(v)
#[1] 1 0 2 0 3 0
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} 
