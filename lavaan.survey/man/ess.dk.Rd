\name{ess.dk}
\alias{ess.dk}
\docType{data}
\title{
Subset of four variables and interviewer number from the 
European Social Survey, round 4 (2008) Danish sample. 
}
\description{
	Consists of indices of political and social trust, political efficacy, and 
	fear of crime. The interviewer number was also recorded. 
	
	The ESS round 4 in Denmark was a face-to-face probability survey 
	with simple random sampling.
}
\usage{data(ess.dk)}
\format{
  A data frame with 1546 observations of 6 variables.
  \describe{
    \item{\code{idno}}{Respondent identifier}
    \item{\code{intnum}}{Interviewer identifier}
    \item{\code{socialTrust}}{Sum of three 0-10 social trust items.}
    \item{\code{efficacy}}{Sum of two 1-5 items measuring "internal political efficacy".}
    \item{\code{systemTrust}}{Sum of three 0-10 political trust items.}
    \item{\code{fearCrime}}{Sum of three 1-4 items measuring fear of crime.}
   }
}
\source{
This dataset was retrieved from \url{http://ess.nsd.uib.no/ess/round4/} and converted to an R dataset.
}
\references{
	Jowell, R., Roberts, C., Fitzgerald, R., & Eva, G. (2007). 
		Measuring attitudes cross-nationally: Lessons from the european social survey. SAGE.
}
\seealso{
\code{\link{lavaan.survey}}
}
\examples{
data(ess.dk)
head(ess.dk)
}
