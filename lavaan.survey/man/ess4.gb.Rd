\name{ess4.gb}
\alias{ess4.gb}
\docType{data}
\title{
  British attitudes towards the welfare state.
}
\description{
European Social Survey (ESS) data from the 2008 (fourth) round in the United Kingdom. The data are from a questionnaire on "what the responsibilities of governments should or should not be". These were factor-analyzed by Roosma, Gelissen, and van Oorschot (2013). Also included are complex survey design variables.
}
\usage{data(ess4.gb)}
\format{
  A data frame with 2273 observations of 13 variables.
  \describe{
    \item{\code{idno}}{Respondent identifier.}
    \item{\code{psu}}{Primary sampling unit (PSU).}
    \item{\code{dweight}}{ESS design weights.}
    \item{\code{stratval}}{Stratification variable (UK regions).}
    
    \item{\code{gvjbevn}}{Job for everyone, governments' responsibility (0-10).}    
    \item{\code{gvhlthc}}{Health care for the sick, governments' responsibility (0-10).}
    \item{\code{gvslvol}}{Standard of living for the old, governments' responsibility (0-10).}
    \item{\code{gvslvue}}{Standard of living for the unemployed, governments' responsibility (0-10).}
    \item{\code{gvcldcr}}{Child care services for working parents, governments' responsibility (0-10).}
    \item{\code{gvpdlwk}}{Paid leave from work to care for sick family, governments' responsibility (0-10).}
    \item{\code{sbprvpv}}{Social benefits/services prevent widespread poverty (1-5).}
    \item{\code{sbeqsoc}}{Social benefits/services lead to a more equal society (1-5).}
    \item{\code{sbcwkfm}}{Social benefits/services make it easier to combine work and family (1-5).}
  }
}
\source{
  This dataset was retrieved from \url{http://www.europeansocialsurvey.org/data/download.html?r=4} and converted to an R dataset.
}
\references{
  Jowell, R., Roberts, C., Fitzgerald, R., & Eva, G. (2007). Measuring attitudes 
    cross-nationally: Lessons from the european social survey. SAGE.

  Oberski, D.L. (2014). lavaan.survey: An R Package for Complex Survey Analysis
    of Structural Equation Models. Journal of Statistical Software, 57(1), 1-27.
    \url{http://www.jstatsoft.org/v57/i01/}.

  Roosma F., Gelissen J., van Oorschot W. (2013). "The Multidimensionality of 
    Welfare State Attitudes: A European Cross-National Study." Social Indicators 
    Research, 113(1), 235-255.
}
\seealso{
  \code{\link{lavaan.survey}}
}
\examples{
  data(ess4.gb)
  
  # Two-factor model based on Roosma et al (2013).
  model.cfa <-    
    "range =~ gvjbevn + gvhlthc + gvslvol + gvslvue + gvcldcr + gvpdlwk
     goals =~ sbprvpv  +  sbeqsoc  +  sbcwkfm"
  
  # Fit the model using lavaan
  fit.cfa.ml <- lavaan(model.cfa, data = ess4.gb, estimator = "MLM",
    meanstructure = TRUE, int.ov.free = TRUE, auto.var = TRUE, 
    auto.fix.first = TRUE, auto.cov.lv.x = TRUE)
  fit.cfa.ml
  
  # Define the complex survey design for ESS 4 in the UK
  des.gb <- svydesign(ids = ~psu, strata = ~stratval, weights = ~dweight, 
    data = ess4.gb)
  
  # Fit the two-factor model while taking the survey design into account.
  fit.cfa.surv <- lavaan.survey(fit.cfa.ml, survey.design = des.gb)
  fit.cfa.surv
}
\author{
  Daniel Oberski - \url{http://daob.nl/} - \email{daniel.oberski@gmail.com}
}
