"stepback" <-

function (y = y, d = dis, alfa = 0.05) 
{
    lm1 <- lm(y ~ ., data = d)
    result <- summary(lm1)
    max <- max(result$coefficients[, 4][-1], na.rm = TRUE)
    while (max > alfa) {
        varout <- names(result$coefficients[, 4])[result$coefficients[, 
            4] == max]
        pos <- position(matrix = d, vari = varout)
        d <- d[, -pos]
        if (length(result$coefficients[, 4][-1]) == 2) {
            min <- min(result$coefficients[, 4][-1], na.rm = TRUE)
            lastname <- names(result$coefficients[, 4])[result$coefficients[, 
                4] == min]
        }
        if (is.null(dim(d))) {
            d <- as.data.frame(d)
            colnames(d) <- lastname
        }
        lm1 <- lm(y ~ ., data = d)
        result <- summary(lm1)
        max <- max(result$coefficients[, 4][-1], na.rm = TRUE)
        if (length(result$coefficients[, 4][-1]) == 1) {
            max <- result$coefficients[, 4][-1]
            if (max > alfa) {
                max = 0
                lm1 <- lm(y ~ 1)
            }
        }
    }
    return(lm1)
}



