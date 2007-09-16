"stepfor" <-
function (y = y, d = d, alfa = 0.05) 
{
    pval <- NULL
    design <- NULL
    j = 1
    resul0 <- summary(lm(y ~ ., data = d))$coefficients[, 4]
    d <- as.data.frame(d[, names(resul0)[-1]])
    for (i in 1:ncol(d)) {
        sub <- cbind(design, d[, i])
        sub <- as.data.frame(sub)
        lm2 <- lm(y ~ ., data = sub)
        result <- summary(lm2)
        pval[i] <- result$coefficients[, 4][j + 1]
    }
    min <- min(pval)
    while (min < alfa) {
        b <- pval == min
        c <- c(1:length(pval))
        pos <- c[b]
        pos <- pos[!is.na(pos)][1]
        design <- cbind(design, d[, pos])
        design <- as.data.frame(design)
        colnames(design)[j] <- colnames(d)[pos]
        j = j + 1
        d <- as.data.frame(d[, -pos])
        pval <- NULL
        if (ncol(d) != 0) {
            for (i in 1:ncol(d)) {
                sub <- cbind(design, d[, i])
                sub <- as.data.frame(sub)
                lm2 <- lm(y ~ ., data = sub)
                result <- summary(lm2)
                pval[i] <- result$coefficients[, 4][j + 1]
            }
            min <- min(pval, na.rm = TRUE)
        }
        else min <- 1
    }
    if (is.null(design)) {
        lm1 <- lm(y ~ 1)
    }
    else {
        lm1 <- lm(y ~ ., data = design)
    }
    return(lm1)
}
