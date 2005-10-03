"PlotProfiles" <-
function(data, cond, main = NULL, cex.xaxis = 0.5, ylim = NULL, repvect, sub = NULL)
{
  pos.vline = rank(repvect, ties.method = "first")
  if (is.null(ylim))   {
    ylim = c(min(as.matrix(data),na.rm = TRUE)*1.1, max(as.matrix(data),na.rm = TRUE)*1.1)
  }
  if(!is.vector(data)) {
    n = dim(data)[2]
    m = dim(data)[1]
    if (m == 1) nom <- rownames(data)
    else  nom <- NULL
    plot(x = c(1:n), y = data[1,], type = "l", col = 1, ylim = ylim, ylab = "expression value", xlab=" ", main = paste("Cluster ", main, "(", m ," genes )",nom), xaxt = "n")
    axis(1,at = 1:n,labels = substr(cond,1,26), cex.axis = cex.xaxis, las = 2)
    abline(v = pos.vline, col = "light gray")
    for (i in 1:dim(data)[1]) {
      lines(x = c(1:n), y = data[i,],col = i)
    }
  } else {
    n = length(data)
    plot(x = c(1:n), y = data, type = "l", col = 1, ylim = ylim, ylab = "expression value", sub ,xaxt = "n", xlab=" ")
    axis(1,at = 1:n, labels = cond, cex.axis = cex.xaxis, las = 2)
    abline(v = pos.vline, col = "light gray")
  }
}

