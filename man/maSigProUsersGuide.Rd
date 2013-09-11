\name{maSigProUsersGuide}
\alias{maSigProUsersGuide}

\title{View maSigPro User's Guide}
\description{
  Finds the location of the maSigPro User's Guide and opens it.
}

\usage{maSigProUsersGuide()}
\arguments{
  \item{view}{logical, to specify if the document is opened using the default PDF document reader.
}

\details{


The function \code{vignette("maSigPro")} will find the short maSigPro Vignette which describes how to obtain the maSigPro User's Guide. 
The User's Guide is not itself a true vignette because it is not automatically generated using Sweave during the package build process. 
This means that it cannot be found using vignette, hence the need for this special function.

If the operating system is other than Windows, then the PDF viewer used is that given by Sys.getenv("R_PDFVIEWER"). 
The PDF viewer can be changed using Sys.putenv(R_PDFVIEWER=). 

}

\examples{
maSigProUsersGuide()
maSigProUsersGuide(view=FALSE)

}
\keyword{UsersGuide}
