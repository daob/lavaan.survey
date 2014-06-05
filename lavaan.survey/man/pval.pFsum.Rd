\name{pval.pFsum}
\alias{pval.pFsum}
\title{F test for model fit of complex survey structural equation models}
\description{
Uses the \code{\link{pFsum}} function from the \pkg{survey} package to obtain a p-value for the overall model fit of a \code{\link{lavaan}} fit object using an F reference distribution, where the denominator degrees of freedom is the design degrees of freedom, \code{\link{degf}}(survey.design).  
}
\usage{
	pval.pFsum(lavaan.fit, survey.design, method = "saddlepoint") 
}
\arguments{
  \item{lavaan.fit}{
	A \code{\linkS4class{lavaan}} object resulting from a lavaan call. 
	
	Since this is the estimator that will be used in the complex sample
	estimates, for comparability it can be convenient to use the same estimator in the call
	generating the \pkg{lavaan} fit object as in the \code{lavaan.survey} call. By default
	this is "MLM".
}
  \item{survey.design}{
	An  \code{\link{svydesign}} object resulting from a call to 
	\code{svydesign} in the \pkg{survey} package. This allows for incorporation of
	clustering, stratification, unequal probability weights, finite
	population correction, and multiple imputation. 
	See the survey documentation for more information.
}
   \item{method}{
		The method by which the distribution of the overall model fit statistic is approximated. See Details on the \code{\link{pFsum}} help page.
   }
}
\details{
With a small number of primary sampling units (design degrees of freedom), the asymptotic chi-square approximation to the distribution of the test statistic may not be entirely accurate. In this case instead of a chi-square, an F-reference distribution using the design degrees of freedom may be used. 

When \code{\link{degf}} is infinite, the p-value output by this function be equal the Satterthwaite ("MLMVS") p-value (see \code{\link{lavaan}} "test" options).

The eigenvalues of the U.Gamma matrix, which is used by \pkg{lavaan} to calculate Satorra-Bentler scaling corrections, will be the coefficients in the 
 mixture of chi-squares distribution (Skinner, Holt & Smith, pp. 86-87).

An anonymous reviewer for the Journal of Statistical Software suggested that 
  "in surveys with small numbers of primary sampling units this sort of 
   correction has often improved the 
   behavior of tests in other contexts."
}
\value{
	A p-value for the overall F test of model fit, adjusted for nonnormality and the complex sampling design.
}
\references{
Skinner C, Holt D, Smith T (1989). Analysis of Complex Surveys. 
  John Wiley & Sons, New York.    

Oberski, D.L. (2014). lavaan.survey: An R Package for Complex Survey Analysis
  of Structural Equation Models. Journal of Statistical Software, 57(1), 1-27.
  \url{http://www.jstatsoft.org/v57/i01/}.
}
\author{
	Daniel Oberski - \url{http://daob.org} - \email{daniel.oberski@gmail.com}
}
\note{Thanks are due to an anonymous reviewer for the Journal of Statistical Software for suggesting this function, and to Yves Rosseel for adjusting the lavaan code to pass along the U.Gamma eigenvalues to the fit object (GitHub commit 225fab0).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{cardinale}}
  \code{\link{lavaan.survey}}
	\code{\link{pFsum}}
	\code{\link{degf}}
	\code{\link{svydesign}}
	\code{\link{lavaan}}
}
\examples{
# Load HolzingerSwineford1939 data
data("HolzingerSwineford1939")

# Create 43 fake clusters
HolzingerSwineford1939$fake.cluster <- rep(1:43, each=7)
# Create survey design object
des <- svydesign(ids=~fake.cluster, probs=~1, data=HolzingerSwineford1939)

# Show the design degrees of freedom: number of clusters - 1
degf(des) # 42

# A reduced factor model that has a larger p-value :
HS.model <- ' visual  =~ x2 + x3
  textual =~ x4 + x5 + x6'

# Fit the factor model without taking complex sampling into account
fit <- cfa(HS.model, data=HolzingerSwineford1939, estimator="MLMVS")

# Fit the factor model, taking the 43 clusters into account
fit.svy <- lavaan.survey(fit, survey.design=des, estimator="MLMVS")

# Calculate the F test p-value. 
# Since degf is only 42, there is a difference with Satterthwaite chi-square
pval.pFsum(fit.svy, survey.design=des)  # 0.0542468133

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{models}
\keyword{regression}
\keyword{robust}
\keyword{multivariate}
