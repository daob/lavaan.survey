\name{lavaan.survey-package}
\alias{lavaan.survey-package}
\docType{package}
\title{
	Complex survey analysis of structural equation models (SEM)
}
\description{
	Takes a lavaan fit object and a complex survey design object as input
	and returns a structural equation modeling analysis based on the fit 
	object, where the complex sampling design is taken into account. 
	
	The structural equation model parameter estimates and standard errors
	are design-based. See Satorra and Muth\'en (1995) for details on the 
	procedure.
}
\details{
\tabular{ll}{
Package: \tab lavaan.survey\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-10-25\cr
License: \tab GPL(>=2)\cr
}

}
\author{
    Author/maintainer: Daniel Oberski <daniel.oberski@gmail.com>
}
\references{
	Oberski, D. and Saris, W. (2012). A model-based procedure to evaluate
	    the relative effects of different TSE components on structural equation
	    model parameter estimates. Presentation given at the International
	    Total Survey Error Workshop in Santpoort, the Netherlands. 
	    \url{http://daob.org/}

	Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics
		and standard errors in covariance structure analysis. 

	Satorra, A., and Muthen, B. (1995). Complex sample data in structural
	   equation modeling. Sociological methodology, 25, 267-316.
}
\keyword{ package }
\seealso{
    \code{\link{lavaan.survey}}
    
    \code{\link{svydesign}}
    \code{\link{svyvar}}

    \code{\link{lavaan}}
    \code{\link{cfa}}
    \code{\link{sem}}
}
\examples{
	# See lavaan.survey function.
}
