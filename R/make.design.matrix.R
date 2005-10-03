"make.design.matrix" <-
function (edesign, degree = 2, time.col = 1,  repl.col = 2,  group.cols = c(3:ncol(edesign))) 
{

## Dummy information
  control.label <- colnames(edesign[,group.cols])[1] 
  if (dim(as.matrix(edesign))[2]>3) {
    dummy.cols <- group.cols[2:length(group.cols)]
    treatm.label <- paste(colnames(edesign)[dummy.cols], "vs", control.label, sep = "")
    groups.label <- c(control.label, treatm.label)
    matrix.dummy <- as.matrix(edesign[,dummy.cols])
    ## Shared time point
    aa <- edesign[,group.cols]
    share <- apply(aa,1,sum)
    if  (!is.element(TRUE, share>1)) {
      dummy <- matrix.dummy
      colnames(dummy) <- treatm.label
    } else  dummy = NULL
## Prepare matrix
    time <- as.matrix(edesign[,time.col])
    colnames(time) <- colnames(edesign)[time.col]
    dis <- cbind(time, dummy)
    rownames(dis) <- rownames(edesign)
    groups.vector <- c(control.label, colnames(dummy))
    colnames.dis <- colnames(dis)
## Lineal terms
    dis <- cbind(dis, dis[,1]*matrix.dummy)
    colnames(dis) <- c(colnames.dis, paste(colnames(edesign)[time.col],"x", colnames(edesign)[dummy.cols], sep=""))
    groups.vector <- c(groups.vector, treatm.label)   
## Generate polynome
    if (degree >= 2){
      for (i in 2 : degree) {
        colnames.dis <- colnames(dis)
        dis <- cbind(dis, edesign[,time.col]^i, edesign[,time.col]^i* edesign[,dummy.cols])
        colnames(dis) <- c(colnames.dis, paste(colnames(edesign)[time.col],i, sep = ""), paste(colnames(edesign)[time.col], "", i, "x", 					colnames(edesign)[dummy.cols],sep = ""))
        groups.vector <- c(groups.vector, groups.label)
      }
    }	
  } else {
## Design matrix when only one group
    dis <- as.matrix(edesign[,time.col])
    colnames(dis) <- colnames(edesign)[time.col]
    rownames(dis) <- rownames(edesign)
    for (i in 2 : degree) {
      colnames.dis <- colnames(dis)
      dis <- cbind(dis, edesign[,time.col]^i) 
      colnames(dis) <- c(colnames.dis, paste(colnames(edesign)[time.col],i, sep = ""))		
    }
    groups.vector <- rep(colnames(edesign)[group.cols], degree)
  }
## Write out results
  output <- list(dis, groups.vector, edesign)
  names(output) <- c("dis", "groups.vector", "edesign")
  output
}