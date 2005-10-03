"p.vector" <-
function (data, design = NULL, Q = 0.05, fdr.method = "BH", min.obs = 3)
{

## Find design matrix
  if (is.data.frame(design) || is.matrix(design)) {
    dis <- design
    groups.vector = NULL
    edesign = NULL
  } else if (is.list(design)) {
    dis <- as.data.frame(design$dis)
    groups.vector <- design$groups.vector
    edesign <- design$edesign
  }
  dat <- as.matrix(data)
  G <- nrow(dat) # initial number of genes
## Min observations
  count.na <- function(x) (length(x)-length(x[is.na(x)]))
  dat <- dat[apply(dat, 1, count.na) >= min.obs,]
## Definitions
  g <- dim(dat)[1] # number of genes to fit
  n <- dim(dat)[2] # number of experimental conditions
  p <- dim(dis)[2] # number of variables
## Initialize variables
  p.vector <- vector(mode = "numeric", length = g)
  p.vector.corre <- vector(mode = "numeric", length = g)
  F.vector <- vector(mode = "numeric", length = g)
## Compute fit per gen
  for (i in 1:g) { 
    dis <-round(dis,5)
    y <- as.numeric(dat[i,])
    reg<-lm(y~.,data = dis)  
    result <- summary(reg)
    if(is.null(result$fstatistic[1])) 
      p.value = 1
    if(!is.null(result$fstatistic[1])) 
      p.vector[i] <- 1-pf(result$fstatistic[1],result$fstatistic[2],result$fstatistic[3])
## Counter
    div <- c(1 : round(g/100)) * 100
    if (is.element (i,div))
      print(paste(c("fitting gene",i, "out of", g), collapse = " "))  
    }
    p.vector <- as.matrix(p.vector)
    rownames(p.vector) <- rownames(dat)
    colnames(p.vector) <- c("p.value")
## Apply FDR
    if (fdr.method == "BH") {
      sortp <- sort(p.vector)
      i <- length(sortp)
      while(i >= 1 && sortp[i] > (Q * i)/g) i <- i - 1
      if(i == 0) {
        alfa <-0
      } else {
       alfa <- sortp[i]
      }
    } else stop ("Only BH fdr.method is currently supported") 
## Compute matrix of selected genes
  SELEC <- dat[which(p.vector <= alfa),] # selection of significant genes
  if (nrow(SELEC)  == 0) print("no significant genes")
## Write out output
  output <-list(SELEC, 
                 p.vector, 
                 G, 
                 g, 
                 alfa, 
                 i, 
                 dis, 
                 dat, 
                 min.obs, 
                 Q, 
                 groups.vector, 
                 edesign)
  names(output) <- c("SELEC", 
                     "p.vector", 
                     "G", 
                     "g", 
                     "alfa", 
                     "i", 
                     "dis", 
                     "dat", 
                     "min.obs", 
                     "Q", 
                     "groups.vector", 
                     "edesign")
  output

}

