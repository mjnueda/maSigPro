"see.genes" <-
function (data, edesign = data$edesign, time.col = 1, repl.col = 2, 
    group.cols = c(3:ncol(edesign)), names.groups = colnames(edesign)[3:ncol(edesign)], 
    cluster.data = 1, groups.vector = data$groups.vector, k = 9, m = 1.45, 
    cluster.method = "hclust", distance = "cor", agglo.method = "ward", 
    show.fit = FALSE, dis = NULL, step.method = "backward", min.obs = 3, 
    alfa = 0.05, nvar.correction = FALSE, show.lines = TRUE, iter.max = 500, 
    summary.mode = "median", color.mode = "rainbow", cexlab = 1, legend = TRUE, 
    newX11 = TRUE,  ylim = NULL, main = NULL, ...) 
{
    time = edesign[, time.col]
    repvect = edesign[, repl.col]
    groups = edesign[, group.cols]
    narrays <- length(time)
    if (!is.null(dim(data))) {
        dat <- as.data.frame(data)
        clusterdata <- data
    }
    else {
        clusterdata <- data[[cluster.data]]
        dat <- as.data.frame(data$sig.profiles)
    }
    clusterdata <- clusterdata
    if (nrow(dat) > 1) {
        dat <- as.data.frame(dat[, (ncol(dat) - length(time) + 
            1):ncol(dat)])
        count.na <- function(x) length(x[is.na(x)])
        NAs <- apply(as.matrix(dat), 1, count.na)
        count.noNa <- function(x) (length(x) - length(x[is.na(x)]))
        dat <- dat[which(apply(as.matrix(dat), 1, count.noNa) >= 
            2), ]
    }
    else {
        NAs <- 1
    }
    kdata <- NULL
    out <- TRUE

    if (nrow(dat) > 1) {
        if (cluster.data != 1 || cluster.data != "sig.profiles") {
            if (any(is.na(clusterdata))) 
                clusterdata[is.na(clusterdata)] <- 0
        }
        else if (is.na(all(dist(clusterdata) > 0)) || (cluster.method == 
            "kmeans" & any(is.na(clusterdata))) || (distance == 
            "cor" & any(sd(t(clusterdata), na.rm = TRUE) == 0))) {
            if (!is.null(kdata)) {
                clusterdata <- kdata
            }
            else {
                clusterdata <- NULL
            }
        }
        clusterdata <- clusterdata
        if (!is.null(clusterdata)) {
            k <- min(k, nrow(dat), na.rm = TRUE)
            if (cluster.method == "hclust") {
                if (distance == "cor") {
                  dcorrel <- matrix(rep(1, nrow(clusterdata)^2), 
                    nrow(clusterdata), nrow(clusterdata)) - cor(t(clusterdata), 
                    use = "pairwise.complete.obs")
                  clust <- hclust(as.dist(dcorrel), method = agglo.method)
                  c.algo.used = paste(cluster.method, "cor", 
                    agglo.method, sep = "_")
                }
                else {
                  clust <- hclust(dist(clusterdata, method = distance), 
                    method = agglo.method)
                  c.algo.used = paste(cluster.method, distance, 
                    agglo.method, sep = "_")
                }
                cut <- cutree(clust, k = k)
            }
            else if (cluster.method == "kmeans") {
                cut <- kmeans(clusterdata, k, iter.max)$cluster
                c.algo.used = paste("kmeans", k, iter.max, sep = "_")
            }
                        else if (cluster.method == "mfuzz") { 
                library(Mfuzz)
                n<-dim(clusterdata)[2]
		    clusterdata[is.na(clusterdata)]<-0
		    temp <- tempfile()
		    write.table(clusterdata, temp, quote = FALSE, sep = "\t", row.names =TRUE, col.names = TRUE)
		    signif <- readExpressionSet(temp)
	 	    cl <- mfuzz(signif, c = k, m = m)
		    clus<-acore(signif,cl=cl,min.acore=(1/k))
		    for(i in 1:k){
		    	clus[[i]]<-transform(clus[[i]],cluster= i )
		    }
		    cut0<-clus[[1]][,c(1,3)]
		    for(i in 2:k){
		    	cut0<-rbind(cut0, clus[[i]][,c(1,3)])
		    }
		    cut<-transform(clusterdata, name="")
		    cut<-transform(cut, cluster=0)
		    cut<-cut[,c(n+1,n+2)]
		    cut[,1]<-rownames(cut)
		    for(i in 1:dim(clusterdata)[1]){
			cut[i,2]<-cut0[cut[i,1],2]
		    }
		    cut<-cut[,2]
                c.algo.used = paste("mfuzz", k, m, sep = "_")
            }
            else stop("Invalid cluster algorithm")
            if (newX11) 
                X11()
            groups <- as.matrix(groups)
            colnames(groups) <- names.groups
            if (k <= 4) 
                par(mfrow = c(2, 2))
            else if (k <= 6) 
                par(mfrow = c(3, 2))
            else if (k > 6) 
                par(mfrow = c(3, 3))
            for (i in 1:(k)) {
                PlotProfiles(data = dat[cut == i, ], repvect = repvect, 
                  main = i, ylim = ylim, color.mode = color.mode, 
                  cond = rownames(edesign), ...)
            }
            if (newX11) 
                X11()
            if (k <= 4) {
                par(mfrow = c(2, 2))
                cexlab = 0.6
            }
            else if (k <= 6) {
                par(mfrow = c(3, 2))
                cexlab = 0.6
            }
            else if (k > 6) {
                par(mfrow = c(3, 3))
                cexlab = 0.35
            }
            for (j in 1:(k)) {
                PlotGroups(data = dat[cut == j, ], show.fit = show.fit, 
                  dis = dis, step.method = step.method, min.obs = min.obs, 
                  alfa = alfa, nvar.correction = nvar.correction, show.lines = show.lines, time = time, 
                  groups = groups, repvect = repvect, summary.mode = summary.mode, 
                  xlab = "time", main = paste("Cluster", j, sep = " "), 
                  ylim = ylim, cexlab = cexlab, legend = legend, 
                  groups.vector = groups.vector, ...)
            }
        }
        else {
            print("warning: impossible to compute hierarchical clustering")
            c.algo.used <- NULL
            cut <- 1
        }
    }
    else if (nrow(dat) == 1) {
        if (newX11) 
            X11()
        PlotProfiles(data = dat, repvect = repvect, main = NULL, 
            ylim = ylim, color.mode = color.mode, cond = rownames(edesign), 
            ...)
        if (newX11) 
            X11()
        PlotGroups(data = dat, show.fit = show.fit, dis = dis, 
            step.method = step.method, min.obs = min.obs, alfa = alfa, nvar.correction = nvar.correction,
            show.lines = show.lines, time = time, groups = groups, 
            repvect = repvect, summary.mode = summary.mode, xlab = "time", 
            main = main, ylim = ylim, cexlab = cexlab, legend = legend, 
            groups.vector = groups.vector, ...)
        c.algo.used <- NULL
        cut <- 1
    }
    else {
        print("warning: NULL data. No visualization possible")
        c.algo.used <- NULL
        cut <- NULL
    }
    OUTPUT <- list(cut, c.algo.used, groups)
    names(OUTPUT) <- c("cut", "cluster.algorithm.used", "groups")
    OUTPUT
}
