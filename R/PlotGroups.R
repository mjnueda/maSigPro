"PlotGroups" <-
function (data,
          edesign = NULL, 
          time = edesign[,1], 
          groups = edesign[,c(3:ncol(edesign))], 
          repvect = edesign[,2], 
          show.fit = FALSE,
          dis = NULL,
          step.method = "backward",
          min.obs = 2,
          alfa = 0.01,
          show.lines = TRUE,
          groups.vector = NULL, 
          xlab = "time", 
          cex.xaxis = 1, 
          ylim = NULL, 
          main = NULL, 
          cexlab = 0.8, 
          legend = TRUE,
          sub = NULL)
{ 
## Compute ratio medio and ylab
  if (!is.vector (data)) {
    yy <- apply (as.matrix(data),2,median, na.rm = TRUE)
    ylab = "averaged expression value"
    if (dim(data)[1] == 1 ) {
      sub <- rownames(data) 
    } else {
      sub <- paste("Mean profile of ", nrow(data), " genes")
    }
  } else if (length(data) != 0 ) {
    yy <- as.numeric(data)
    sub <-  rownames(data)
    ylab = "expression value"
  } else stop ("empty data")
## Calculate groups
  if (is.null(ncol(groups))) {
    ncol = 1
    legend = FALSE
    codeg = "group"
  } else  {
    ncol = ncol(groups)
    codeg <- as.character(colnames(groups))
}

##  Average over time replicates	
  reps <- i.rank(repvect)
  y <- vector(mode = "numeric", length = length(unique(reps)))
  x <- vector(mode = "numeric", length = length(unique(reps)))
  g <- matrix(nrow = length(unique(reps)), ncol = ncol)
  for (k in 1:length(y)) {
    y[k] <- mean(yy[reps == k], na.rm = TRUE)
    x[k] <- mean(time[reps == k])
    for(j in 1:ncol) {
      g[k,j] <-mean(groups[reps == k,j])
    }
  }  
## Plot
    if (is.null(ylim))  ylim = c(min(as.numeric(yy),na.rm = TRUE), max(as.numeric(yy),na.rm = TRUE))
    abcissa <- x
    xlim = c(min(abcissa,na.rm = TRUE), max(abcissa,na.rm = TRUE)*1.3)
    color1 <- as.numeric(sort(factor(colnames(groups)))) + 1
    color2 <- groups
    for (j in 1:ncol) {
      color2[,j] <- color2[,j] * j } 
    color2 <- as.vector(apply(color2, 1, sum) + 1)
#### Plot data points
    plot(x = time, 
         y = yy, 
         pch = 21, 
         xlab = xlab, 
         ylab = ylab, 
         xaxt = "n", 
         main = main, 
         sub = sub, 
         ylim = ylim, 
         xlim = xlim, 
         cex = cexlab,
	 col= color2)
    axis(1,at = unique(abcissa),labels = unique(abcissa), cex.axis = cex.xaxis)
    if(show.fit) {
      rm <- matrix(yy, nrow = 1, ncol = length(yy))
      rownames(rm) <- c("ratio medio")
      colnames(rm) <- rownames(dis)
      fit.y <- T.fit(rm, design = dis, step.method = step.method,  min.obs = min.obs, alfa = alfa)
      betas <- fit.y$coefficients
  }
#### Compute and plot regression curves
  for (i in 1: ncol(groups)) {
    group <- g[,i]
    if ((show.fit) && !is.null(betas)) {
      li <- c(2:6)
      a <- reg.coeffs (coefficients = betas, groups.vector = groups.vector, group = colnames(groups)[i])
      a <- c(a, rep(0,(5-length(a))))
      curve(a[1] + a[2]*x + a[3]*(x^2) + a[4]*(x^3) + a[5]*(x^4), from = xlim[1], to = xlim[2], col = color1[i],  
            add = TRUE, lty = li[i])
    }
#### Plot data curves
    if (show.lines) {
      lx <- abcissa[group != 0]
      ly <- y[group != 0]
      ord <- order(lx)
      lxo <- lx[ord]
      lyo <- ly[ord]
      lines(lxo, lyo, col = color1[i])      
    }
  }
  op <- par(bg = "white")   
  if (legend) 
    legend(max(abcissa,na.rm = TRUE) * 1.02, ylim[1],legend = codeg, text.col = color1, col = color1, cex = cexlab, lty = 1, yjust = 0)
  par(op)
}

