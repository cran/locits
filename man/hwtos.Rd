\name{hwtos}
\alias{hwtos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Haar wavelet test for (second-order) stationarity for arbitrary length time series.

%%  ~~function to do ... ~~
}
\description{NOTE: CURRENTLY THIS FUNCTION IS NOT INCLUDED IN THE PACKAGE.
	USE hwtos2. This function computes the raw wavelet periodogram of the
	arbitrary time series vector \code{x}. The periodogram is then
	subject to a hypothesis test to see if its expectation over time,
	for different scales, is constant. The constancy test is carried out
	using tests on its Haar wavelet coefficients. The overall test is
	for second-order stationarity (e.g. constant variance, constant
	acf function, mean is assumed zero).
}
\usage{
hwtos(x, alpha = 0.05, lowlev = 1, WTscale = NULL, maxSD = NULL,
    verbose = FALSE, silent = FALSE, UseCForVarip2 = TRUE, OPLENGTH = 1e+05,
    mc.method = p.adjust.methods)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The time series you wish to test for second-order stationarity.
	Minimum length series that this function will operate for is 20.
	However, for short series the power of the test might not be good
	and could be investigated via simulation that reflect your particular
	circumstances.
	This should be a stochastic series. The function
        will report an error if \code{x} is a constant function. The
        function might not work properly if it contains a significant trend
        or patches of non-stochastic observations.
}
  \item{alpha}{The (nominal) size of the hypothesis test.
}
  \item{lowlev}{Controls the lowest scale of the wavelet periodogram that
	gets analyzed. Generally, leave this parameter alone.
}
  \item{WTscale}{Controls the finest scale of the Haar wavelet transform
	of a particular wavelet periodogram scale. Generally, we have to
	stay away from the finest Haar wavelet transform scales of the
	periodogram as the test relies on a central limit theorem effect
	which only "kicks in" when the Haar wavelet scale is medium-to-coarse.
	Generally, leave this argument alone.
}
  \item{maxSD}{Parameter which controls which scales go towards overall
	variance calculation. Generally, leave alone.
}
  \item{verbose}{If \code{TRUE} then informative messages are printed.
	If \code{FALSE} only limited informational messages are printed
	unless \code{silent=TRUE}.
}
  \item{silent}{If \code{TRUE} then no messages are printed during the
	operation of the function at all.
}
  \item{UseCForVarip2}{If \code{TRUE} then fast C code is used for the
	variance calculation, otherwise slower R code is used. 
}
  \item{OPLENGTH}{Some of the internal functions require workspace to
	perform their calculations. In exceptional circumstances more
	static workspace might be required and so this argument might
	need to be higher than the default. However, the code will tell
	you how high this number will need to be. The code can, with
	default arguments, handle series that are up to 30000 in length.
	However, at 35000 the OPLENGTH parameter will need to be increased.
}
  \item{mc.method}{Method to control overall size for test taking into
	account multiple comparisons. The default argument is
	\code{p.adjust.methods} which is the same as the default argument
	to the \code{p.adjust} function in R. This includes a number
	of the popular methods such as "Holm", "Bonferonni" and "FDR",
	for example.
	}
}
\details{This function computes all possible Haar wavelet coefficients
	of the time series \code{x}. Then, squares those to obtain the
	raw wavelet periodogram. Then the test of stationarity works
	by taking each level of the raw wavelet periodogram and subjecting
	it to another (decimated) Haar wavelet transform and then assessing
	whether any of those coefficients is significantly different to
	zero. It does this by using a Gaussian approximation first
	introduced by Neumann and von Sachs (2000). This is a multiple
	testing problem: many individual wavelet coefficients need to
	be assessed simultaneously and the user can choose the type of
	assessment using the \code{mc.method} argument.
}
\value{An object of class \code{tosANYN}. This is a list containing the
	following components.
\item{nreject}{The number of wavelet coefficients that reject the null
	hypothesis of being zero.}
\item{mc.method}{The multiple comparison method used.}
\item{AllTS}{All the t-statistics. This is a list containing J levels,
	where J is the number of periodogram levels. Each slot in
	the \code{AllTS} list itself contains a Haar wavelet transform
	object (\code{hwtANYN}) which are the t-statistics associated with
	each Haar wavelet coefficient of the Haar raw wavelet periodogram.
	}
\item{AllPVal}{As \code{AllTS} but for p-values}
\item{alpha}{The size of the test}
\item{x}{The time series that was analyzed}
\item{xSD}{The estimated mean spectrum value for each level of the
	spectrum, mean over time that is.}
\item{allTS}{A vector containing all of the test statistics. So,
	the information in \code{AllTS} but arranged as a single vector}
\item{allpvals}{As \code{allTS} but for p-values. These values have been
	adjusted to take account of the multiple comparisons. See the
	vaector \code{allpvals.unadjust} for an unadjusted set.}
\item{allbigscale}{The wavelet periodogram scale associated with each
	t-statistic in \code{allTS}.}
\item{alllitscale}{As for \code{allbigscale} but for the wavelet transform
	of the wavelet periodogram.}
\item{allindex}{As for \code{allbigscale} but the wavelet coefficient index
	in the Haar wavelet transform of the wavelet periodogram}
\item{alllv}{The maximum number of wavelet coefficients in a particular
	Haar wavelet scale of a particular scale of the wavelet periodogram.
	Note, this information is useful because the wavelet transforms
	are computed on arbitrary length objects and so keeping track of
	the number of coefficients per scale is useful later, e.g. for
	plotting purposes. This information is not required in the
	dyadic case because the coefficient vector lengths are completely
	predictable.}
\item{allpvals.unadjust}{A vector of p-values that has not been adjusted
	by a multiple hypothesis test technique.} 
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


\seealso{\code{link{hwt}},
	\code{\link{hwtos2}},
	\code{\link{plot.tosANYN}},
	\code{\link{print.tosANYN}},
	\code{\link{summary.tosANYN}}
}
\examples{
#
# Generate test data set of non-dyadic length
#
v3 <- rnorm(300)
#
# Run the test of stationarity
#
\dontrun{v3.TOS <- hwtos(v3)}
#
#Scales get printed 
#8    7    6    5    4    3    2 
#
\dontrun{print(v3.TOS)}
#Class 'tosANYN' : Stationarity Object for Arbitrary Length Data :
#	~~~~~~~  : List with 14 components with names
#	nreject mc.method AllTS AllPVal alpha x xSD allTS
#	allpvals allbigscale alllitscale allindex alllv
#	allpvals.unadjust 
#
#
#summary(.):
#----------
#There are  54  hypothesis tests altogether
#There were  0  reject(s)
#P-val adjustment method was:  holm 
#
# Note, nothing got rejected. So accept the H_0 null hypothesis of stationarity.
# This is precisely what you'd expect operating on iid Gaussians.
#
# Let's construct obvious example of non-stationarity.
#
v4 <- c(rnorm(150), rnorm(150,sd=3))
#
# I.e. v4 has sharp variance change halfway along
# Now compute test of stationarity
#
\dontrun{v4.TOS <- hwtos(v4)}
#
#  Print out results
#
\dontrun{print(v4.TOS)}
#
#Class 'tosANYN' : Stationarity Object for Arbitrary Length Data :
#       ~~~~~~~  : List with 14 components with names
#       nreject mc.method AllTS AllPVal alpha x xSD allTS
#       allpvals allbigscale alllitscale allindex alllv
#       allpvals.unadjust 
#
#
#summary(.):
#----------
#There are  54  hypothesis tests altogether
#There were  5  reject(s)
#P-val adjustment method was:  holm 
#Listing rejects...
#P:  7  HWTlev:  2  Max Poss Ix:  2  Indices:  2 
#P:  7  HWTlev:  1  Max Poss Ix:  1  Indices:  1 
#P:  6  HWTlev:  1  Max Poss Ix:  1  Indices:  1 
#P:  5  HWTlev:  1  Max Poss Ix:  1  Indices:  1 
#P:  4  HWTlev:  1  Max Poss Ix:  1  Indices:  1 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
