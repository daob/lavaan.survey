\name{pisa.be.2003}
\alias{pisa.be.2003}
\docType{data}
\title{
  Belgian school childrens' math efficacy and measured ability in the PISA study.
}
\description{
  Data from the OECD's 2003 Programme for International Student Assessment (PISA) in Belgium. 

 Also included are 80 so-called "replicate weights" calculated by Westat which must be accounted for using complex survey procedures.
}
\usage{data(pisa.be.2003)}
\format{
  A data frame with 8796 observations of 102 variables.
  \describe{
    \item{\code{PV1MATH1}}{"Plausible values" for child's overall math ability (imputation 1).}
    \item{\code{PV1MATH2}}{"Plausible values" for child's overall math ability (imputation 2).}
    \item{\code{PV1MATH3}}{"Plausible values" for child's overall math ability (imputation 3).}
    \item{\code{PV1MATH4}}{"Plausible values" for child's overall math ability (imputation 4).}
    \item{\code{ST31Q01}}{Feel confident doing task: "timetable" 1 (very) - 4 (not at all).}
    \item{\code{ST31Q02}}{Feel confident doing task: "discount" (1-4).}
    \item{\code{ST31Q03}}{Feel confident doing task: "area" (1-4).}
    \item{\code{ST31Q04}}{Feel confident doing task: "graphs" (1-4).}
    \item{\code{ST31Q05}}{Feel confident doing task: "linear" (1-4).}
    \item{\code{ST31Q06}}{Feel confident doing task: "distance" (1-4).}
    \item{\code{ST31Q07}}{Feel confident doing task: "quadratics" (1-4).}
    \item{\code{ST31Q08}}{Feel confident doing task: "rate" (1-4).}
    
    \item{\code{ST32Q02}}{"I am just not good at Mathematics" 1 (strongly agree) - 4 (strongly disagree).}
    \item{\code{ST32Q04}}{"I get good marks in Mathematics" (1-4).}
    \item{\code{ST32Q06}}{"I learn Mathematics quickly" (1-4).}
    \item{\code{ST32Q07}}{"I have always believed that Mathematics is one of my best subjects" (1-4).}
    \item{\code{ST32Q09}}{"In my Mathematics class, I understand even the most difficult work" (1-4).}
    \item{\code{ESCS}}{Index of Socio-Economic and Cultural Status.}
    \item{\code{male}}{Gender (1=Female, 2=male).}
    \item{\code{school.type}}{Difficulty level of secondary studies (1-3).}
    \item{\code{W_FSTUWT}}{Overall survey weight for students.}
    \item{\code{W_FSTR1}}{BRR replicate.}
\item{\code{W_FSTR2}}{BRR replicate.}
\item{\code{W_FSTR3}}{BRR replicate.}
\item{\code{W_FSTR4}}{BRR replicate.}
\item{\code{W_FSTR5}}{BRR replicate.}
\item{\code{W_FSTR6}}{BRR replicate.}
\item{\code{W_FSTR7}}{BRR replicate.}
\item{\code{W_FSTR8}}{BRR replicate.}
\item{\code{W_FSTR9}}{BRR replicate.}
\item{\code{W_FSTR10}}{BRR replicate.}
\item{\code{W_FSTR11}}{BRR replicate.}
\item{\code{W_FSTR12}}{BRR replicate.}
\item{\code{W_FSTR13}}{BRR replicate.}
\item{\code{W_FSTR14}}{BRR replicate.}
\item{\code{W_FSTR15}}{BRR replicate.}
\item{\code{W_FSTR16}}{BRR replicate.}
\item{\code{W_FSTR17}}{BRR replicate.}
\item{\code{W_FSTR18}}{BRR replicate.}
\item{\code{W_FSTR19}}{BRR replicate.}
\item{\code{W_FSTR20}}{BRR replicate.}
\item{\code{W_FSTR21}}{BRR replicate.}
\item{\code{W_FSTR22}}{BRR replicate.}
\item{\code{W_FSTR23}}{BRR replicate.}
\item{\code{W_FSTR24}}{BRR replicate.}
\item{\code{W_FSTR25}}{BRR replicate.}
\item{\code{W_FSTR26}}{BRR replicate.}
\item{\code{W_FSTR27}}{BRR replicate.}
\item{\code{W_FSTR28}}{BRR replicate.}
\item{\code{W_FSTR29}}{BRR replicate.}
\item{\code{W_FSTR30}}{BRR replicate.}
\item{\code{W_FSTR31}}{BRR replicate.}
\item{\code{W_FSTR32}}{BRR replicate.}
\item{\code{W_FSTR33}}{BRR replicate.}
\item{\code{W_FSTR34}}{BRR replicate.}
\item{\code{W_FSTR35}}{BRR replicate.}
\item{\code{W_FSTR36}}{BRR replicate.}
\item{\code{W_FSTR37}}{BRR replicate.}
\item{\code{W_FSTR38}}{BRR replicate.}
\item{\code{W_FSTR39}}{BRR replicate.}
\item{\code{W_FSTR40}}{BRR replicate.}
\item{\code{W_FSTR41}}{BRR replicate.}
\item{\code{W_FSTR42}}{BRR replicate.}
\item{\code{W_FSTR43}}{BRR replicate.}
\item{\code{W_FSTR44}}{BRR replicate.}
\item{\code{W_FSTR45}}{BRR replicate.}
\item{\code{W_FSTR46}}{BRR replicate.}
\item{\code{W_FSTR47}}{BRR replicate.}
\item{\code{W_FSTR48}}{BRR replicate.}
\item{\code{W_FSTR49}}{BRR replicate.}
\item{\code{W_FSTR50}}{BRR replicate.}
\item{\code{W_FSTR51}}{BRR replicate.}
\item{\code{W_FSTR52}}{BRR replicate.}
\item{\code{W_FSTR53}}{BRR replicate.}
\item{\code{W_FSTR54}}{BRR replicate.}
\item{\code{W_FSTR55}}{BRR replicate.}
\item{\code{W_FSTR56}}{BRR replicate.}
\item{\code{W_FSTR57}}{BRR replicate.}
\item{\code{W_FSTR58}}{BRR replicate.}
\item{\code{W_FSTR59}}{BRR replicate.}
\item{\code{W_FSTR60}}{BRR replicate.}
\item{\code{W_FSTR61}}{BRR replicate.}
\item{\code{W_FSTR62}}{BRR replicate.}
\item{\code{W_FSTR63}}{BRR replicate.}
\item{\code{W_FSTR64}}{BRR replicate.}
\item{\code{W_FSTR65}}{BRR replicate.}
\item{\code{W_FSTR66}}{BRR replicate.}
\item{\code{W_FSTR67}}{BRR replicate.}
\item{\code{W_FSTR68}}{BRR replicate.}
\item{\code{W_FSTR69}}{BRR replicate.}
\item{\code{W_FSTR70}}{BRR replicate.}
\item{\code{W_FSTR71}}{BRR replicate.}
\item{\code{W_FSTR72}}{BRR replicate.}
\item{\code{W_FSTR73}}{BRR replicate.}
\item{\code{W_FSTR74}}{BRR replicate.}
\item{\code{W_FSTR75}}{BRR replicate.}
\item{\code{W_FSTR76}}{BRR replicate.}
\item{\code{W_FSTR77}}{BRR replicate.}
\item{\code{W_FSTR78}}{BRR replicate.}
\item{\code{W_FSTR79}}{BRR replicate.}
\item{\code{W_FSTR80}}{BRR replicate.}
    \item{\code{WVARSTRR}}{Randomized final variance stratum (1-80).}
}
}
\source{This dataset was retrieved from \url{http://www.oecd.org/pisa/} and converted to an R dataset.

For more information regarding the variables and the model estimated below, please see Ferla et al. (2009) and/or the PISA manual (OECD, 2009). 
}
\references{
Ferla J, Valcke M, Cai Y (2009). "Academic Self-Efficacy and Academic Self-
  Concept: Reconsidering Structural Relationships." Learning and Individual 
  Differences, 19(4), 499-505.

Oberski, D.L. (2014). lavaan.survey: An R Package for Complex Survey Analysis
  of Structural Equation Models. Journal of Statistical Software, 57(1), 1-27.
  \url{http://www.jstatsoft.org/v57/i01/}.

OECD (2009). PISA Data Analysis Manual: SPSS and SAS. 2nd edition. OECD.
}
\seealso{
  \code{\link{lavaan.survey}}
}
\examples{
  data(pisa.be.2003)


  # Simplified version of Ferla et al. (2009) model.
  model.pisa <- "
      math =~ PV1MATH1 + PV1MATH2 + PV1MATH3 + PV1MATH4
      neg.efficacy =~ ST31Q01 + ST31Q02 + ST31Q03 + ST31Q04 + 
                      ST31Q05 + ST31Q06 + ST31Q07 + ST31Q08
      neg.selfconcept =~ ST32Q02 + ST32Q04 + ST32Q06 + ST32Q07 + ST32Q09
  
      neg.selfconcept ~ neg.efficacy + ESCS + male
      neg.efficacy ~ neg.selfconcept + school.type + ESCS + male
      math ~ neg.selfconcept + neg.efficacy + school.type + ESCS + male
  "
  
  # Fit the model using lavaan
  fit <- lavaan(model.pisa, data = pisa.be.2003, auto.var = TRUE, std.lv = TRUE,
    meanstructure = TRUE, int.ov.free = TRUE, estimator = "MLM")
  
  ## Not run due to CRAN policies, 
  ##    uncomment below to account for replicate weights:
  
  ## Define the survey design using the BRR replicate weights provided by PISA
  ##   Note that these settings will work for _any_ analysis of PISA data...
  #des.rep <- svrepdesign(ids = ~1, weights = ~W_FSTUWT, data = pisa.be.2003, 
  #  repweights = "W_FSTR[0-9]+", type = "Fay", rho = 0.5)
  
  ## Fit the SEM model accounting for replicate weights
  #fit.surv <- lavaan.survey(lavaan.fit = fit, survey.design =  des.rep)
  
  #fit # Show fitmeasures results
  #fit.surv # Show fitmeasures results
}
\author{
  Daniel Oberski - \url{http://daob.nl/} - \email{daniel.oberski@gmail.com}
}
