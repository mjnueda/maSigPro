"maSigPro" <-
function (data, 
          edesign,  
          matrix ="AUTO", 
          groups.vector = NULL, 
          degree = 2, 
          time.col = 1, 
          repl.col = 2,   
          group.cols = c(3:ncol(edesign)), 
          Q = 0.05,
          alfa = Q,
          step.method = "backward", 
          rsq = 0.7,
          min.obs = 3,
          vars = "groups",
          cluster.data = 1, 
          add.IDs = FALSE,  
          IDs = NULL, 
          matchID.col = 1, 
          only.names = FALSE,
          k = 9,
          cluster.method = "kmeans", 
          distance = "cor", 
          agglo.method = "complete", 
          iter.max = 500, 
          trat.repl.spots = "none",  
          index = IDs[,(matchID.col+1)],
          match= IDs[,matchID.col], 
          rs=0.7, 
          show.fit = TRUE,
          show.lines = TRUE, 
          pdf = TRUE,
          cexlab = 0.8,
          legend = TRUE,
          main = NULL, ...)
{

## Calculate design matrix
  if (matrix == "AUTO") {
    print ("running design")
    design <- make.design.matrix(edesign = edesign, 
                                 degree = degree, 
                                 time.col = time.col, 
                                 repl.col =  repl.col, 
                                 group.cols = group.cols)
    dis <- design$dis
    groups.vector <- design$groups.vector
  } else {
    design = dis <- matrix
    groups.vector <- groups.vector
  }
  cluster.algorithm <- NULL 
  groups <- unique(groups.vector)
  STOP = FALSE
##  Compute global fit
  print ("running p.vector")
  fit <- p.vector(data = data, design = design, Q = Q , min.obs = min.obs)
  if (is.null(fit$SELEC) || nrow(fit$SELEC) == 0 )  {
    summary <- c("no significant genes")
    print("maSigPro halted at p.vector")
    output<- list(summary,
                  fit$dat,
                  fit$G, 
                  edesign, 
                  dis, 
                  fit$min.obs,
                  fit$p.vector, 
                  Q)
    names(output) <- c("summary",
                       "input.data", 
                       "G", 
                       "edesign", 
                       "dis", 
                       "min.obs",
                       "p.vector", 
                       "Q")
    STOP = TRUE
  }

  if(!STOP) {
## Compute stepwise fit
    print ("running T.fit")   
    tstep <- T.fit (data = fit , step.method = step.method , min.obs = min.obs, alfa = alfa)
    if (is.null(tstep$sol) || nrow(tstep$sol) == 0 ) {
      summary <- c("no significant genes")
      print("maSigPro halted at tstep")
      output<- list(summary,
                    fit$dat,
                    fit$G, 
                    edesign, 
                    dis, 
                    fit$min.obs,
                    fit$p.vector, 
                    tstep$variables, 
                    tstep$g, 
                    fit$alfa, 
                    step.method, 
                    Q, 
                    alfa,
                    tstep$influ.info)
      names(output) <- c("summary",
                         "input.data", 
                         "G", 
                         "edesign", 
                         "dis", 
                         "min.obs",
                         "p.vector", 
                         "variables", 
                         "g", 
                         "p.vector.alfa",
                         "step.method", 
                         "Q", 
                         "step.alfa",
                         "influ.info")					
      STOP=TRUE
    }
    if(!STOP) {
## Get sig.genes
      print ("running get.siggenes")
      got.genes <- get.siggenes (tstep, 
                                 vars = vars, 
                                 rsq = rsq,
                                 groups.vector = groups.vector, 
                                 add.IDs = add.IDs, 
                                 IDs = IDs, 
                                 matchID.col = matchID.col, 
                                 only.names = only.names, 
                                 trat.repl.spots = trat.repl.spots, 
                                 index = index, 
                                 match = match)
      summary <- got.genes$summary
      sig.genes <- got.genes$sig.genes
      sig.genes <-sig.genes	
## Graphical display
      if (!is.null(sig.genes)) {
        if (pdf) {
          if  (!is.null(main)) {
            pdf (file = paste(main, "pdf", sep = "."), title = main)
          } else {
            pdf (file = "Results.pdf")
          }	
        }
      if (!only.names) {
        if(vars != "all") {
          for(i in 1:length(sig.genes)) {
            if (nrow(sig.genes[[i]][[1]]) > 0) {
              print (paste ("running see.genes ", i))
              cluster <- see.genes(data = sig.genes[[i]], 
                                   cluster.data = cluster.data, 
                                   k = k, 
                                   cluster.method = cluster.method, 
                                   distance = distance, 
                                   agglo.method = agglo.method,
                                   show.fit = show.fit, 
                                   dis = dis,
                                   step.method = step.method,
                                   min.obs = min.obs,
                                   alfa = alfa,
                                   show.lines = show.lines,
                                   cexlab = cexlab,
                                   newX11 = FALSE,
                                   legend = legend,
                                   main = paste(main, names(sig.genes[i]), sep=" "), ...)
              sig.genes[[i]][[1]] <- cbind(sig.genes[[i]][[1]], cluster$cut)
              cluster.algorithm <- cluster$cluster.algorithm.used
              groups <- cluster$groups
            } 
          }
        } else {
          if (nrow(sig.genes[[1]]) > 0) {
            print ("running see.genes")
            cluster <- see.genes(data = sig.genes,
                                 cluster.data = cluster.data, 
                                 k = k, 
                                 cluster.method = cluster.method, 
                                 distance = distance , 
                                 agglo.method = agglo.method,
                                 show.fit = show.fit,
                                 dis = dis,
                                 step.method = step.method,
                                 min.obs = min.obs,
                                 alfa = alfa,
                                 show.lines = show.lines,
                                 cexlab = cexlab,
                                 legend = legend,
                                 newX11 = FALSE,
                                 main = main, ...)
            sig.genes[[1]] <- cbind(sig.genes[[1]], cluster$cut)
            cluster.algorithm <- cluster$cluster.algorithm.used
            groups <- cluster$groups
          }
        }            
      }
      dev.off()
    } else print("maSigPro halted at get.siggenes")

## Write out results
    output<- list(summary,	
                  sig.genes, 
                  fit$dat,
                  fit$G, 
                  edesign, 
                  dis, 
                  fit$min.obs,
                  fit$p.vector, 
                  tstep$variables, 
                  tstep$g, 
                  fit$alfa, 
                  step.method, 
                  Q, 
                  alfa,
                  tstep$influ.info,
                  vars, 
                  cluster.algorithm, 
                  groups)
    names(output) <- c("summary",
                       "sig.genes", 
                       "input.data", 
                       "G", 
                       "edesign", 
                       "dis", 
                       "min.obs",
                       "p.vector", 
                       "variables", 
                       "g", 
                       "p.vector.alfa",
                       "step.method", 
                       "Q", 
                       "step.alfa",
                       "influ.info", 
                       "select.vars", 
                       "cluster.algorithm.used", 
                       "groups")
    }
  }
  output
}

