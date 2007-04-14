"suma2Venn" <-
function (x, ...) 
{

library("limma")
total <- NULL
## Number of different elements
  for (i in 1:ncol(x)) 
    total <- unique(c(total, as.character(x[,i])))
total <- total[total!=" "]
## Creat matrix of counts
  counts <- matrix(0,nrow = length(total), ncol = ncol(x))
  colnames(counts) <- colnames(x)
  rownames(counts) <- total
  for(j in 1:ncol(x)) {
        counts[, j] [is.element(total, x[,j])] <- 1
  }
## create Venn diagram
  if (ncol(counts) <= 3) {
    vennDiagram(counts, ...)
  }  else print("venn diagram cannot be made: too many columns")
  
}

