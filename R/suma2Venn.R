"suma2Venn" <-
function (x, ...) 
{

library("limma")
total <- NULL
## Number of different elements
  for (i in 1:ncol(x)) 
    total <- c(total, unique(as.character(x[,i])))
## Creat matrix of counts
  counts <- matrix(0,nrow = length(total), ncol = ncol(x))
  colnames(counts) <- colnames(x)
  for(j in 1:ncol(x)) {
    for (i in 1 : length(total)){
      if (is.element(total[i], x[,(j)]))
        counts[i,j] <- 1
    }
  }
  counts
## create Venn diagram
  if (ncol(counts) >= 3) {
    vennDiagram(counts, ...)
  }  else print("venn diagram cannot be made: too many columns")
  
}

