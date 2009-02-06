\name{edesign.OD}
\alias{edesign.OD}
\docType{data}
\title{ Experimental design with a measured independent variable}
\description{
  \code{edesign.OD} contains the experimental design of a E.coli growth time course microarray experiment with a temperature shift treatment. The OD of each culture was measured and used in the experimental design as independent variable.
}
\usage{data(edesign.OD)}
\format{
  A data frame with 52 rows and the following 4 variables.
  \describe{
    \item{\code{OD}}{a numeric vector. Indicates the OD value of the sampled culture}
    \item{\code{Replicate}}{a numeric vector}
    \item{\code{37}}{a numeric vector. No temperature shitf treatment}
    \item{\code{SHIFT}}{a numeric vector. Temperature shift treatment}
  }
}

\examples{
data(edesign.OD)
## maybe str(edesign.OD) ; plot(edesign.OD) ...
}
\keyword{datasets}
