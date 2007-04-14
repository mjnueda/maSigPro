\name{maSigPro}
\alias{maSigPro}

\title{Wrapping function for identifying significant differential gene expression profiles in micorarray time course experiments }
\description{
     \code{maSigPro} performs a whole maSigPro analysis for a times series gene expression experiment.
      The function sucesively calls the functions \code{\link{make.design.matrix}}(optional), \code{\link{p.vector}}, \code{\link{T.fit}},
      \code{\link{get.siggenes}} and \code{\link{see.genes}}.
}
\usage{
maSigPro(data, edesign, matrix = "AUTO", groups.vector = NULL, 
    degree = 2, time.col = 1, repl.col = 2, group.cols = c(3:ncol(edesign)), 
    Q = 0.05, alfa = Q, nvar.correction = FALSE, step.method = "backward", rsq = 0.7,
    min.obs = 3, vars = "groups", significant.intercept = "dummy", cluster.data = 1, 
    add.IDs = FALSE, IDs = NULL, matchID.col = 1, only.names = FALSE, k = 9, 
    cluster.method = "hclust", distance = "cor", agglo.method = "ward", iter.max = 500, 
    summary.mode = "median", color.mode = "rainbow", trat.repl.spots = "none",
    index = IDs[, (matchID.col + 1)], match = IDs[, matchID.col], rs = 0.7, 
    show.fit = TRUE, show.lines = TRUE, pdf = TRUE, cexlab = 0.8, 
    legend = TRUE, main = NULL, ...)
}

\arguments{
  \item{data}{ matrix with normalized gene expression data. Genes must be in rows and arrays in columns. Row names must contain geneIDs
    \newline \quad	(argument of \code{\link{p.vector}})}
  \item{edesign}{ matrix of experimental design. Row names must contain arrayIDs
    \newline \quad      (argument of \code{\link{make.design.matrix}} and \code{\link{see.genes}})}
  \item{matrix}{ design matrix for regression analysis. By default design is calculated with make.design.matrix
     \newline \quad     (argument of \code{\link{p.vector}} and \code{\link{T.fit}}, by default computed by \code{\link{make.design.matrix}})}
  \item{groups.vector}{ vector indicating experimental group of each variable 
    \newline \quad      (argument of \code{\link{get.siggenes}} and \code{\link{see.genes}}, by default computed by \code{\link{make.design.matrix}})}
  \item{degree}{ the degree of the regression fit polynome. \code{degree} = 1 returns lineal regression, \code{degree} = 2 returns quadratic regression, etc...
    \newline \quad      (argument of \code{\link{make.design.matrix}}) }
  \item{time.col}{ column in edesign containing time values. Default is first column
    \newline \quad      (argument of \code{\link{make.design.matrix}} and \code{\link{see.genes}}) }
  \item{repl.col}{ column in edesign containing coding for replicates arrays. Default is second column
    \newline \quad      (argument of \code{\link{make.design.matrix}} and \code{\link{see.genes}}) }
  \item{group.cols}{columns in \code{edesign} indicating the coding for each group of the experiment (see \code{\link{make.design.matrix}})
    \newline \quad      (argument of \code{\link{make.design.matrix}} and \code{\link{see.genes}}) }
  \item{Q}{ level of false discovery rate (FDR) control
    \newline \quad      (argument of \code{\link{p.vector}})}
  \item{alfa}{ significance level used for variable selection in the stepwise regression
    \newline \quad      (argument of \code{\link{T.fit}})}
  \item{nvar.correction}{logical for indicating correcting of stepwise regression significance level
    \newline \quad      (argument of \code{\link{T.fit}})}
  \item{step.method}{ argument to be passed to the step function. 
    \newline \quad  Can be either \code{"backward"}, \code{"forward"}, \code{"two.ways.backward"} or \code{"two.ways.forward"} }
  \item{rsq}{ cut-off level at the R-squared value for the stepwise regression fit. 
    \newline \quad Only genes with R-squared greater than \code{rsq} are selected }
  \item{min.obs}{ genes with less than this number of true numerical values will be excluded from the analysis 
    \newline \quad      (argument of \code{\link{p.vector}} and \code{\link{T.fit}})}
  \item{vars}{ variables for which to extract significant genes
    \newline \quad     (argument of \code{\link{get.siggenes}})}
  \item{significant.intercept}{ experimental groups for which significant intercept coefficients are considered
    \newline \quad     (argument of \code{\link{get.siggenes}})}
  \item{cluster.data}{ Type of data used by the cluster algorithm
    \newline \quad     (argument of \code{\link{see.genes}}) }
  \item{add.IDs}{ logical indicating whether to include additional gene id's in the significant genes result 
    \newline \quad     (argument of \code{\link{get.siggenes}}) }
  \item{IDs}{ matrix contaning additional gene id information (required when \code{add.IDs} is TRUE)
    \newline \quad     (argument of \code{\link{get.siggenes}}) }
  \item{matchID.col}{ number of matching column in matrix IDs for adding genes ids
    \newline \quad     (argument of\code{\link{get.siggenes}}) }
  \item{only.names}{ logical. If TRUE, expression values are ommited in the significant genes result
    \newline \quad     (argument of \code{\link{get.siggenes}}) }
  \item{k}{ number of clusters 
    \newline \quad(argument of \code{\link{see.genes}}) }
  \item{cluster.method}{ clustering method for data partioning
    \newline \quad     (argument of \code{\link{see.genes}}) }
  \item{distance}{ distance measurement function used when \code{cluster.method} is \code{"hclust"}
    \newline \quad     (argument of \code{\link{see.genes}}) }
  \item{agglo.method}{ aggregation method used when \code{cluster.method} is \code{"hclust"}
    \newline \quad     (argument of \code{see.genes}) }
  \item{iter.max}{ number of iterations when \code{cluster.method} is \code{"kmeans"}
    \newline \quad     (argument of \code{\link{see.genes}}) }
  \item{summary.mode}{ the method to condensate expression information when more than one gene is present in the data. 
    \newline \quad Possible values are \code{"representative"} and \code{"median"} 
    \newline \quad     (argument of \code{\link{PlotGroups}}) }  
 \item{color.mode}{ color scale for plotting profiles. Can be either \code{"rainblow"} or \code{"gray"}
    \newline \quad     (argument of \code{\link{PlotProfiles}}) }  
  \item{trat.repl.spots}{ treatment givent to replicate spots. Possible values are \code{"none"} and \code{"average"}
    \newline \quad     (argument of \code{\link{get.siggenes}}) }
  \item{index}{ argument of the \code{\link{average.rows}} function to use when \code{trat.repl.spots} is \code{"average"}
    \newline \quad     (argument of \code{\link{get.siggenes}}) }
  \item{match}{ argument of the \code{link{\average.rows}} function to use when \code{trat.repl.spots} is \code{"average"}
    \newline \quad     (argument of \code{\link{get.siggenes}}) }
  \item{rs}{ minimun pearson correlation coefficient for replicated spots profiles to be averaged
    \newline \quad     (argument of \code{\link{get.siggenes}}) }
  \item{show.fit}{ logical indicating whether regression fit curves must be plotted
    \newline \quad     (argument of \code{\link{see.genes}}) }
  \item{show.lines}{ logical indicating whether a line must be drawn joining plotted data points for reach group
    \newline \quad     (argument of \code{\link{see.genes}})}
  \item{pdf}{ logical indicating whether a pdf results file must be generated
    \newline \quad     (argument of \code{\link{see.genes}})}
  \item{cexlab}{ graphical parameter maginfication to be used for x labels in plotting functions }
  \item{legend}{ logical indicating whether legend must be added when plotting profiles 
   \newline \quad (argument of \code{\link{see.genes}}) }
  \item{main}{ title for pdf results file }
  \item{...}{ other graphical function arguments}
}
\details{

    maSigPro finds and display genes with significant profile differences in time series gene expression experiments.
    The main, compulsory, input parameters for this function are a matrix of gene expression data (see \code{\link{p.vector}} for details)
    and a matrix describing  experimental design (see \code{\link{make.design.matrix}} or \code{\link{p.vector}} for details). In case extended
    gene ID information is wanted to be included in the result of significant genes, a third IDs matrix containing this 
    information will be required (see \code{\link{get.siggenes}} for details).

    Basiscally in the function calls subsequent steps of the maSigPro approach which is:
  \newline \qquad Make a general regression model with dummies to indicate different experimental groups.
  \newline \qquad Select significant genes on the basis of this general model, applying fdr control.
  \newline \qquad Find significant variables for each gene, using stepwise regression.
  \newline \qquad Extract and display significant genes for any set of variables or experimental groups.
}
\value{
  \item{summary}{a vector or matrix listing significant genes for the variables given by the function parameters}
  \item{sig.genes}{a list with detailed information on the significant genes found for the variables given by the function parameters. Each element of the list is also a list containing:
  \newline \qquad \code{sig.profiles}: expression values of significant genes.The cluster assingment of each gene is given in the last column
  \newline \qquad \code{coefficients}: regression coefficients for significant genes
  \newline \qquad \code{t.score}: value of the t statistics of significant genes
  \newline \qquad \code{sig.pvalues}: p-values of the regression coefficients for significant genes
  \newline \qquad \code{g}: number of genes
  \newline \qquad \dots :arguments passed by previous functions}
  \item{input.data }{input analysis data}
  \item{G}{ number of input genes}
  \item{edesign}{matrix of experimental design}
  \item{dis}{regression design matrix}
  \item{min.obs}{imputed value for minimal number of true observations}
  \item{p.vector}{vector containing the computed p-values of the general regression model for each gene}
  \item{variables }{variables in the general regression model}
  \item{g}{number of signifant genes}
  \item{p.vector.alfa }{p-vlaue at FDR = \code{Q} control}
  \item{step.method}{imputed step method for stepwise regression}
  \item{Q}{imputed value for false discovery rate (FDR) control}
  \item{step.alfa}{inputed significance level in stepwise regression}
  \item{influ.info }{data frame of genes containing influencial data}
}
\references{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talón, T. 2005.
maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments.
 }
\author{ Ana Conesa, aconesa@ivia.es; María José Nueda, mj.nueda@ua.es}
\seealso{ \code{\link{make.design.matrix}}, \code{\link{p.vector}}, \code{\link{T.fit}}, \code{\link{get.siggenes}}, \code{\link{see.genes}}}
\examples{

#### GENERATE TIME COURSE DATA
## generate n random gene expression profiles of a data set with 
## one control plus 3 treatments, 3 time points and r replicates per time point.

tc.GENE <- function(n, r,
             var11 = 0.01, var12 = 0.01,var13 = 0.01,
             var21 = 0.01, var22 = 0.01, var23 =0.01,
             var31 = 0.01, var32 = 0.01, var33 = 0.01,
             var41 = 0.01, var42 = 0.01, var43 = 0.01,
             a1 = 0, a2 = 0, a3 = 0, a4 = 0,
             b1 = 0, b2 = 0, b3 = 0, b4 = 0,
             c1 = 0, c2 = 0, c3 = 0, c4 = 0)
{

  tc.dat <- NULL
  for (i in 1:n) {
    Ctl <- c(rnorm(r, a1, var11), rnorm(r, b1, var12), rnorm(r, c1, var13))  # Ctl group
    Tr1 <- c(rnorm(r, a2, var21), rnorm(r, b2, var22), rnorm(r, c2, var23))  # Tr1 group
    Tr2 <- c(rnorm(r, a3, var31), rnorm(r, b3, var32), rnorm(r, c3, var33))  # Tr2 group
    Tr3 <- c(rnorm(r, a4, var41), rnorm(r, b4, var42), rnorm(r, c4, var43))  # Tr3 group
    gene <- c(Ctl, Tr1, Tr2, Tr3)
    tc.dat <- rbind(tc.dat, gene)
  }
  tc.dat
}

## Create 270 flat profiles
flat <- tc.GENE(n = 270, r = 3)
## Create 10 genes with profile differences between Ctl and Tr1 groups
twodiff <- tc.GENE (n = 10, r = 3, b2 = 0.5, c2 = 1.3)
## Create 10 genes with profile differences between Ctl, Tr2, and Tr3 groups
threediff <- tc.GENE(n = 10, r = 3, b3 = 0.8, c3 = -1, a4 = -0.1, b4 = -0.8, c4 = -1.2)
## Create 10 genes with profile differences between Ctl and Tr2 and different variance
vardiff <- tc.GENE(n = 10, r = 3, a3 = 0.7, b3 = 1, c3 = 1.2, var32 = 0.03, var33 = 0.03)
## Create dataset
tc.DATA <- rbind(flat, twodiff, threediff, vardiff)
rownames(tc.DATA) <- paste("feature", c(1:300), sep = "")
colnames(tc.DATA) <- paste("Array", c(1:36), sep = "")
tc.DATA[sample(c(1:(300*36)), 300)] <- NA  # introduce missing values

#### CREATE EXPERIMENTAL DESIGN
Time <- rep(c(rep(c(1:3), each = 3)), 4)
Replicates <- rep(c(1:12), each = 3)
Control <- c(rep(1, 9), rep(0, 27))
Treat1 <- c(rep(0, 9), rep(1, 9), rep(0, 18))
Treat2 <- c(rep(0, 18), rep(1, 9), rep(0,9))
Treat3 <- c(rep(0, 27), rep(1, 9))
edesign <- cbind(Time, Replicates, Control, Treat1, Treat2, Treat3)
rownames(edesign) <- paste("Array", c(1:36), sep = "")

#### RUN maSigPro
tc.test <- maSigPro (tc.DATA, edesign, degree = 2, vars = "groups", main = "Test")

tc.test$g  # gives number of total significant genes
tc.test$summary  # shows significant genes by experimental groups
tc.test$sig.genes$Treat1$sig.pvalues  # shows pvalues of the significant coefficients 
                                      # in the regression models of the significant genes 
                                      # for Control.vs.Treat1 comparison

}
\keyword{ manip }
\keyword{ models }
