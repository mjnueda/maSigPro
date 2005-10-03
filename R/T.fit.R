"T.fit" <-
function (data, design = data$dis, step.method = "backward",  min.obs = data$min.obs, alfa = data$Q)
{
## Find design matrix
  if (is.list(data)) {
    dat <- as.matrix(data$SELEC)
    dat <- rbind(c(rep(1,ncol(dat))), dat)
    groups.vector <- data$groups.vector
    edesign <- data$edesign
    G <- data$g
  } else {
    G <- nrow(data)
    data <- rbind(c(rep(1,ncol(data))), data)
    dat <- as.matrix(data)
    count.na <- function(x) (length(x) - length(x[is.na(x)]))
    dat <- dat[apply(dat, 1, count.na) >= min.obs,]
    groups.vector = NULL
    edesign = NULL
  }
  dis <<- as.data.frame(design)
## Extract dimensions
  g <- (dim(dat)[1]-1) 
  n <- dim(dat)[2] 
  p <- dim(dis)[2] 
  vars.in <- colnames(dis) 
## Inicialize variables
  sol <- NULL
  coefficients <- NULL
  t.score <- NULL
  sig.profiles <- NULL
  influ.info<-matrix(NA, nrow=nrow(dis), ncol=1)
  rownames(influ.info) <- rownames(dis)
## Compute gene per gene STEP regression
  for (i in 2:(g+1)) { 
    y <- as.numeric(dat[i,])
    name <- rownames(dat)[i]	 
    if (step.method == "backward") {
      reg <- stepback(y = y, d = dis, alfa = alfa)
    } else if (step.method == "forward") {
      reg <- stepfor(y = y, d = dis, alfa = alfa)
    } else if (step.method == "two.ways.backward") {
      reg <- two.ways.stepback(y = y, d = dis, alfa = alfa)
    } else if (step.method == "two.ways.forward") {
      reg <- two.ways.stepfor(y = y, d = dis, alfa = alfa)
    } else stop ("stepwise method must be one of backward, forward, two.ways.backward, two.ways.forward")
## Counter
    div <- c(1 : round(g/100)) * 100
    if (is.element (i,div))
      print(paste(c("fitting gene",i, "out of", g), collapse = " "))
## Find only possible variables
    lmf <- lm(y~., data = as.data.frame(dis))
    result <- summary(lmf)
    novar <-vars.in[!is.element(vars.in,names(result$coefficients[,4]))]     
## Find influence measurements
    influ <- influence.measures(reg)$is.inf
    influ <- cbind(influ[,ncol(influ)-3],influ[,ncol(influ)-1])
    influ1 <- which(apply(influ,1,all))  
    if(length(influ1) != 0) { 
      match <- match(rownames(dis),rownames(influ))
      colnames(influ) <- rep(rownames(dat)[i],2)
      influ.info<-cbind(influ.info,influ[match,])
    }
    result <- summary(reg)
## Extract statistics per gene
    if (!is.null(result$fstatistic[1])) {
      k <- i
      p.value <- 1-pf(result$fstatistic[1], result$fstatistic[2], result$fstatistic[3]) 
      beta.coeff <- result$coefficients[,1]   
      beta.p.valor <- result$coefficients[,4]                                                     
      coeff <- rep(0, (length(vars.in)+1))
      if (length(novar) != 0) {
        for (m in 1:length(novar)) {
          coeff [position(dis,novar[m])+1] <-NA
        }
      }
      p.valor = t <- rep("NA", (length(vars.in)+1))
      vars.out <- NULL
      if (result$coefficients[,4][rownames(result$coefficients) == "(Intercept)"] < alfa) {
        coeff[1] <- result$coefficients[,1][rownames(result$coefficients) == "(Intercept)"]
        p.valor[1] <- result$coefficients[,4][rownames(result$coefficients) == "(Intercept)"]
        t[1] <- result$coefficients[,3][rownames(result$coefficients) == "(Intercept)"]
      }
      for (j in 2 :length(coeff)) {
        if(is.element(vars.in[j-1],rownames(result$coefficients))) {
          coeff[j] <- result$coefficients[,1][rownames(result$coefficients) == vars.in[j-1]]
          p.valor[j] <- result$coefficients[,4][rownames(result$coefficients) == vars.in[j-1]]
          t[j] <- result$coefficients[,3][rownames(result$coefficients) == vars.in[j-1]]
          vars.out <- cbind(vars.out, vars.out[j-1])
        }
      }
      if (!all(is.na(p.valor))) {
        sol <- rbind(sol, as.numeric(c(p.value, result$r.squared, p.valor)))
        coefficients <- rbind(coefficients, coeff)
        t.score <- rbind(t.score, t)
        sig.profiles <- rbind(sig.profiles,y)
        h <- nrow(sol)
        rownames(sol)[h] <- name
        rownames(coefficients)[h] <- name
        rownames(t.score)[h] <- name
        rownames(sig.profiles)[h] <- name
      }                         
    }     
  }
##  Write matrix with results
  if (!is.null(sol)) {
    sol <- as.data.frame(sol)
    coefficients <- as.data.frame(coefficients)
    t.score <- as.data.frame(t.score)
    sig.profiles <- as.data.frame(sig.profiles)
    colnames(sol) <-c ("p-value", "R-squared", "p.valor_beta0",paste("p.valor_",vars.in, sep = "")) 
    colnames(coefficients) <-c ("beta0", paste("beta", vars.in, sep = "")) 
    colnames(t.score) <-c ("t.score_beta0", paste("t.score_", vars.in, sep = "")) 
    colnames(sig.profiles) <- colnames(dat)
  } 	
  influ.info <- influ.info[,-1]
## Write out output
  output <- list(sol, 
		 sig.profiles, 
		 coefficients, 
		 t.score,   
		 vars.in, 
		 G,
		 g, 
		 dat, 
		 dis, 
		 step.method,
		 groups.vector, 
		 edesign, 
		 influ.info)
  names(output) <- c("sol", 
		     "sig.profiles", 
		     "coefficients", 
		     "t.score", 
		     "variables",  
		     "G", 
		     "g", 
		     "dat", 
		     "dis", 
		     "step.method", 
		     "groups.vector", 
		     "edesign",
 		     "influ.info")
  output
}

