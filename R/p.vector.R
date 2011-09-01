"p.vector" <-
function (data, design = NULL, Q = 0.05, MT.adjust = "BH", min.obs = 3) 
{
    if (is.data.frame(design) || is.matrix(design)) {
        dis <- design
        groups.vector = NULL
        edesign = NULL
    }
    else if (is.list(design)) {
        dis <- as.data.frame(design$dis)
        groups.vector <- design$groups.vector
        edesign <- design$edesign
    }
    dat <- as.matrix(data)
    dat <- dat[, as.character(rownames(dis))]
    G <- nrow(dat)
    count.na <- function(x) (length(x) - length(x[is.na(x)]))
    dat <- dat[apply(dat, 1, count.na) >= min.obs, ]
    g <- dim(dat)[1]
    n <- dim(dat)[2]
    p <- dim(dis)[2]
    p.vector <- vector(mode = "numeric", length = g)
    p.vector.corre <- vector(mode = "numeric", length = g)
    F.vector <- vector(mode = "numeric", length = g)
    for (i in 1:g) {
        dis <- round(dis, 5)
        y <- as.numeric(dat[i, ])
        reg <- lm(y ~ ., data = dis)
        result <- summary(reg)
        if (is.null(result$fstatistic[1])) 
            p.value = 1
        if (!is.null(result$fstatistic[1])) 
            p.vector[i] <- 1 - pf(result$fstatistic[1], result$fstatistic[2], 
                result$fstatistic[3])
        div <- c(1:round(g/100)) * 100
        if (is.element(i, div)) 
            print(paste(c("fitting gene", i, "out of", g), collapse = " "))
    }
    p.vector <- as.matrix(p.vector)
    rownames(p.vector) <- rownames(dat)
    colnames(p.vector) <- c("p.value")
    p.adjusted <- p.adjust(p.vector, method = MT.adjust, n = length(p.vector))
    BH.alfa=NULL
    if (MT.adjust == "BH") {
        sortp <- sort(p.vector)
        i <- length(sortp)
        while (i >= 1 && sortp[i] > (Q * i)/g) i <- i - 1
        if (i == 0) {
            BH.alfa <- 0
        }
        else {
            BH.alfa <- sortp[i]
        }
    }

   genes.selected <- rownames(dat)[which(p.adjusted <= Q)]
   SELEC <- as.matrix(as.data.frame(dat)[genes.selected, ])

     if (nrow(SELEC) == 0) 
        print("no significant genes")
    output <- list(SELEC, p.vector, p.adjusted, G, g, BH.alfa, nrow(SELEC), dis, dat, 
        min.obs, Q, groups.vector, edesign)
    names(output) <- c("SELEC", "p.vector", "p.adjusted", "G", "g", "BH.alfa", 
        "i", "dis", "dat", "min.obs", "Q", "groups.vector", "edesign")
    output
}