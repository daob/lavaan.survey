\name{cardinale}
\alias{cardinale}
\docType{data}
\title{Nutrient experiment on patches of algae in Californian streams.}
\description{
An experiment on patches of algae in Californian streams by Cardinale, Bennett, Nelson, and Gross (2009). The authors chose 20 streams in the Mono Lake and Owens River basins in the Sierra Nevada. In each stream, they placed 5 or 10 PVC elbows containing different levels of nutrients and a small patch of agar on which algae could grow. They then returned to the streams about 42 days later and measured 1) species diversity in the stream, 2) species diversity in each patch, 3) biomass of the algae, and 4) rate of oxygen production on each patch. 
}
\usage{data(cardinale)}
\format{
  A data frame with 127 observations of 7 variables.
  
  \describe{
    \item{\code{PatchDiversity}}{Number of species on agar.}
    \item{\code{Biomass}}{Standing algal biomass.}
    \item{\code{O2Production}}{Gross primary oxygen production.}
    \item{\code{logNutrient}}{Experimentally manipulated nutrient supply rate (log_10 transformed).}
    \item{\code{logNutrient2}}{log(nutrient, base=10)^2.}
    \item{\code{StreamDiversity}}{Number of species in stream.}
    \item{\code{Stream}}{The stream in which the experiment was performed.}
  }
}
\source{
 The results of the fit may be compared with Figure 5 in Cardinale et al. (2009, p. 1237).

 Model and data from this example were obtained from Jarrett Byrnes' GitHub:
    https://github.com/jebyrnes/Ecological-SEMs-in-lavaan
    
 Note that I changed the variable names to be easier to understand, e.g.
    logNutrient was logN, PatchDiversity was SA, etc.
}
\references{
Cardinale BJ, Bennett DM, Nelson CE, Gross K (2009). "Does Productivity Drive
   Diversity or Vice Versa? A Test of the Multivariate Productivity-Diversity 
   Hypothesis in Streams." Ecology, 90(5), 1227-1241.

Oberski, D.L. (2014). lavaan.survey: An R Package for Complex Survey Analysis
  of Structural Equation Models. Journal of Statistical Software, 57(1), 1-27.
  \url{http://www.jstatsoft.org/v57/i01/}.
}
\seealso{
  \code{\link{lavaan.survey}}
}
\examples{
  data(cardinale)
  summary(cardinale)
  
  model.card <- '
     PatchDiversity ~ logNutrient + logNutrient2 + StreamDiversity
     Biomass ~ PatchDiversity + logNutrient
     O2Production ~ logNutrient + Biomass
     logNutrient ~~ logNutrient2'
  
  fit.card <- sem(model.card, data = cardinale, fixed.x = FALSE, estimator = "MLM")
  
  des.card <- svydesign(ids = ~Stream, probs = ~1, data = cardinale)
  fit.card.survey <- lavaan.survey(fit.card, des.card, estimator = "MLM")
  
  pval.pFsum(fit.card.survey, survey.design = des.card)

}
\author{
  Daniel Oberski - \url{http://daob.nl/} - \email{daniel.oberski@gmail.com}
}
