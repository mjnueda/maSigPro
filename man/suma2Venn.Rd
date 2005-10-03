\name{suma2Venn}
\alias{suma2Venn}

\title{ Creates a Venn Diagram from a matrix of characters}
\description{
\code{suma2Venn} transforms a matrix of characters into a binary matrix and creates 
a vennDiagram of the common elements between columns
}

\usage{
suma2Venn(x, ...)
}

\arguments{
  \item{x}{ data frame of character values}
  \item{\dots}{ plotting arguments for the vennDiagram function }
}
\details{
  \code{suma2Venn} creates a list of all elements of a matrix or data frame of characters and computes  
the presence/absence of each element in each column of the matrix. This results is a numeric
matrix of 1 and 0 which can be taken by the \code{\link{vennDiagram}} to generate a Venn Plot
}

\value{
\code{suma2Venn} returns a Venn Plot such as that created by the \code{\link{vennDiagram}} funcion
}

\author{ Ana Conesa, aconesa@ivia.es}

\seealso{  \code{\link{vennDiagram}}}
\examples{
a <- c("a","b","c", "d", "e", NA, NA)
b <- c("a","b","f", NA, NA, NA, NA)
c <- c("b","e","f", "h", "i", "j", "k")
x <- cbind(a, b,c)
suma2Venn(x)

}
\keyword{ misc }
\keyword{ aplot }
