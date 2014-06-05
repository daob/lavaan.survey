\name{lavaan.survey}
\alias{lavaan.survey}
\title{
	Complex survey analysis of structural equation models (SEM)
}
\description{
	Takes a lavaan fit object and a complex survey design object as input
	and returns a structural equation modeling analysis based on the fit 
	object, where the complex sampling design is taken into account. 
	
	The structural equation model parameter estimates are "aggregated" (Skinner, Holt & Smith 1989), i.e. they consistently estimate parameters aggregated over any 
  clusters and strata and no explicit modeling of the effects of clusters and strata
  is involved. Standard errors are design-based. 
	See Satorra and Muthen (1995) and references below for details on the procedure.
  
 Both the pseudo-maximum likelihood (PML) procedure popular in the SEM world 
 (e.g. Asparouhov 2005; Stapleton 2006) and 
 weighted least squares procedures similar to aggregate regression modeling with
 complex sampling (e.g. Fuller 2009, chapter 6) are implemented.
 
 It is possible to give a list of multiply imputed datasets to svydesign as data.
 \code{lavaan.survey} will then apply the standard Rubin (1987) formula to obtain
 point and variance estimates under multiple imputation. Some care is required with
 this procedure when survey weights are also involved, however (see Notes).
}
\usage{
lavaan.survey(lavaan.fit, survey.design, 
	     estimator=c("MLM", "MLMV", "MLMVS", "WLS", "DWLS", "ML"),
	     estimator.gamma=c("default","Yuan-Bentler"))
	 

}
%- maybe also 'usage' for other objects documented here.
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
   \item{estimator}{
   		The estimator used determines how parameter estimates are obtained, 
   		how standard errors are calculated, and how the test statistic and 
   		all measures derived from it are adjusted. See \code{\link{lavaan}}.
   		
 		The default estimator is MLM. It is recommended to use one
		of the ML estimators.
   }
   \item{estimator.gamma}{
	Whether to use the usual estimator of Gamma as given by \code{svyvar} (the variance-covariance
	matrix of the observed variances and covariances), or apply some kind
	of smoothing or adjustment. Currently the only other option is the
	Yuan-Bentler (1998) adjustment based on model residuals.
   }
}
\details{
	The user specifies a complex sampling design with the \pkg{survey} package's
	\code{\link{svydesign}} function, and a structural equation model with
	\code{\link{lavaan}}.
	
	 \code{lavaan.survey} follows these steps:
	\enumerate{
		\item The covariance matrix of the observed variables
		(or matrices in the case of multiple
		group analysis) is estimated using the \code{svyvar} command from the
		\pkg{survey} package. 
	       	\item The asymptotic covariance matrix of the variances and
		covariances is obtained from the \code{svyvar} output  (the "Gamma"
		matrix)
		\item The last step depends on the estimation method chosen:
		
		\enumerate{
		\item[MLM, MLMV, MLMVS] The \pkg{lavaan} model is re-fit using Maximum Likelihood
		with the covariance matrix as data. After normal-theory ML
		estimation, the standard errors (\code{vcov} matrix), likelihood ratio
		("chi-square") statistic, and all derived fit indices and
		statistics are adjusted for the complex sampling design using
		the Gamma matrix. I.e. the Satorra-Bentler (SB) corrections are
		obtained ("MLM" estimation in \pkg{lavaan} terminology). This procedure
    is equivalent to "pseudo"-maximum likelihood (PML).
		\item[WLS, DWLS] The \pkg{lavaan} model is re-fit using Weighted Least Squares
		with the covariance matrix as data, and the Moore-Penrose inverse
		of the Gamma matrix as estimation weights. If DWLS is chosen
		only the diagonal of the weight matrix is used.
		}
	}


}
\value{
	An object of class \code{\linkS4class{lavaan}}, where the estimates, 
	standard errors, \code{vcov} matrix, chi-square statistic, and fit measures 
	based on the chi-square take into account the complex survey 
	design. Several methods are available for \code{\linkS4class{lavaan}} 
	objects, including a \code{summary} method.}
\references{
Asparouhov T (2005). Sampling Weights in Latent Variable Modeling. Structural 
  Equation Modeling, 12(3), 411-434.

Bollen, K, Tueller, S, Oberski, DL (2013). Issues in the Structural Equation 
  Modeling of Complex Survey Data. In: Proceedings of the 59th World Statistics 
  Congress 2013 (International Statistical Institute, ed.), Hong Kong.
  \url{http://daob.nl/publications/}

Fuller WA (2009). Sampling Statistics. John Wiley & Sons, New York.

Kim J, Brick J, Fuller WA, Kalton G (2006). On the Bias of the Multiple-Imputation 
  Variance Estimator in Survey Sampling. Journal of the Royal Statistical Society: 
  Series B (Statistical Methodology), 68(3), 509-521.

Oberski, D.L. (2014). lavaan.survey: An R Package for Complex Survey Analysis
  of Structural Equation Models. Journal of Statistical Software, 57(1), 1-27.
  \url{http://www.jstatsoft.org/v57/i01/}.

Oberski, D. and Saris, W. (2012). A model-based procedure to evaluate
    the relative effects of different TSE components on structural equation
    model parameter estimates. Presentation given at the International
    Total Survey Error Workshop in Santpoort, the Netherlands. 
    \url{http://daob.nl/publications/}

Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics
	and standard errors in covariance structure analysis. 

Satorra, A., and Muthen, B. (1995). Complex sample data in structural
   equation modeling. Sociological methodology, 25, 267-316.

Skinner C, Holt D, Smith T (1989). Analysis of Complex Surveys. 
  John Wiley & Sons, New York.

Stapleton L (2006). An Assessment of Practical Solutions for Structural Equation 
  Modeling with Complex Sample Data. Structural Equation Modeling, 13(1), 28-58.

Stapleton L (2008). Variance Estimation Using Replication Methods in Structural 
  Equation Modeling with Complex Sample Data. Structural Equation Modeling, 15(2), 183-210.

Yuan K, Bentler P (1998). Normal Theory Based Test Statistics in Structural Equation
  Modelling. British Journal of Mathematical and Statistical Psychology, 51(2), 289-309.  
}
\author{
	Daniel Oberski - \url{http://daob.nl/} - \email{daniel.oberski@gmail.com}
}
\note{    
    1) Some care should be taken when applying multiple imputation with survey
    weights. The weights should be incorporated in the imputation, and even 
    then the variance produced by the usual Rubin (1987) estimator may not
    be consistent (Kott 1995; Kim et al. 2006).

    If multiple imputation is used to deal with unit nonresponse,
    calibration and/or propensity score weighting with jackknifing may be a 
    more appropriate method. See the \pkg{survey} package.

 2) Note that when using PML or WLS, the Gamma matrix need not be positive definite.
  Preliminary investigations suggest that it often is not. This may happen due to 
  reduction of effective sample size from clustering, for instance. 
  In itself this need not be a problem, depending on the restrictiveness of the model.
  In such cases \code{lavaan.survey} checks explicitly whether the covariance matrix
  of the parameter estimates is still positive definite and produces a warning otherwise. 

    3) Currently only structural equation models for continuous variables are 
    implemented.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{pval.pFsum}}

  \code{\link{cardinale}}
  \code{\link{ess.dk}}
  \code{\link{ess4.gb}}
  \code{\link{liss}}
  \code{\link{pisa.be.2003}}

  \code{\link{svydesign}}
  \code{\link{lavaan}}
}
\examples{
###### A single group example #######

# European Social Survey Denmark data (SRS)
data(ess.dk)

# A saturated model with reciprocal effects from Saris & Gallhofer
dk.model <- "
  socialTrust ~ 1 + systemTrust + fearCrime
  systemTrust ~ 1 + socialTrust + efficacy
  socialTrust ~~ systemTrust
"
lavaan.fit <- lavaan(dk.model, data=ess.dk, auto.var=TRUE, estimator="MLM")
summary(lavaan.fit)

# Create a survey design object with interviewer clustering
survey.design <- svydesign(ids=~intnum, prob=~1, data=ess.dk)

survey.fit <- lavaan.survey(lavaan.fit=lavaan.fit, survey.design=survey.design)
summary(survey.fit)



###### A multiple group example #######

data(HolzingerSwineford1939)

# The Holzinger and Swineford (1939) example - some model with complex restrictions
HS.model <- ' visual  =~ x1 + x2 + c(lam31, lam31)*x3
              textual =~ x4 + x5 + c(lam62, lam62)*x6
              speed   =~ x7 + x8 + c(lam93, lam93)*x9 
             speed ~ textual 
             textual ~ visual'

# Fit multiple group per school
fit <- lavaan(HS.model, data=HolzingerSwineford1939,
              int.ov.free=TRUE, meanstructure=TRUE,
              auto.var=TRUE, auto.fix.first=TRUE, group="school",
              auto.cov.lv.x=TRUE, estimator="MLM")
summary(fit, fit.measures=TRUE)

# Create fictional clusters in the HS data
set.seed(20121025)
HolzingerSwineford1939$clus <- sample(1:100, size=nrow(HolzingerSwineford1939), replace=TRUE)
survey.design <- svydesign(ids=~clus, prob=~1, data=HolzingerSwineford1939)

summary(fit.survey <- lavaan.survey(fit, survey.design))


# For more examples, please see the Journal of Statistical Software Paper, 
#  the accompanying datasets ?cardinale ?ess4.gb ?liss ?pisa.be.2003
#  and my homepage http://daob.nl/ 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{models}
\keyword{regression}
\keyword{robust}
\keyword{multivariate}
