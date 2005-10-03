"get.siggenes" <-
function (tstep, 
          rsq = 0.7, 
          add.IDs = FALSE, 
          IDs = NULL, 
          matchID.col = 1, 
          only.names = FALSE, 
          vars = c("all","each","groups"), 
          groups.vector = NULL, 
          trat.repl.spots = "none", 
          index = IDs[,(matchID.col+1)], 
          match = IDs[,matchID.col],
          r = 0.7) 
{

## Find variables 
  dis <- tstep$dis
  edesign <- tstep$edesign
  groups.vector <- tstep$groups.vector
  not.all.empty <- function(x) (is.element(FALSE,x == " "))
##  Filter genes for the given rsq threshold
  if (any(tstep$sol[,2] > rsq)) {
    sig.pvalues <- tstep$sol[which(tstep$sol[,2] > rsq),] 
    sig.profiles <- tstep$sig.profiles[which(tstep$sol[,2] > rsq),]
    coefficients <- tstep$coefficients[which(tstep$sol[,2] > rsq),]
    t.score <- tstep$t.score[which(tstep$sol[,2] > rsq),]
## Vars == all
      if (vars == "all") { 
        sigs <- sig.profiles
        summary <-  rownames(sig.profiles)                               
        if(only.names) {
          sigs <- rownames(sig.profiles)
        }
        if(add.IDs) { 
          row.names <- rownames(sig.profiles)
          ids <- IDs[is.element(IDs[,matchID.col], rownames(sig.profiles)),]
          if(!only.names) {
            sigs <- cbind(ids, sigs)
          } else {
            sigs <- ids 
          }
          rownames(sigs) <- row.names
        }
        coeffs <-coefficients
        ts <- t.score
        ps <- sig.pvalues
        sig.genes <- list(sigs, coeffs, ts, ps,nrow(sig.profiles), edesign, groups.vector)
        names(sig.genes) <- c("sig.profiles", "coefficients","t.score", 
              "sig.pvalues","g", "edesign", "groups.vector")     		 	
## Vars == each
      } else  if (vars == "each") {           
        sig.genes <- as.list(paste ("var", c("independ", colnames(dis)), sep=".") )
        summary <- matrix(" ", ncol = ncol(dis) + 1, nrow = nrow(sig.profiles))
        colnames (summary) <- c("independ", colnames(dis))
        for (i in 1:ncol(summary)) {
          sigs <- sig.profiles[which(!is.na(sig.pvalues[,(2+i)])),]
          coeffs <- coefficients[which(!is.na(sig.pvalues[,(2+i)])),]
          ts <- t.score[which(!is.na(sig.pvalues[,(2+i)])),]
          ps <- sig.pvalues[which(!is.na(sig.pvalues[,(2+i)])),]
          if (nrow(sigs) > 0) 
            names.sigs <- rownames(sigs) 
          else names.sigs <- NULL
          summary[,i] <- c(names.sigs, rep(" ", nrow(sig.profiles)-nrow(sigs)))
          sig.genes [[i]] <- list(sigs,coeffs, ts, ps,nrow(sigs), edesign, groups.vector)
          names(sig.genes[[i]]) <- c("sig.profiles", "coefficients","t.score", 
                "sig.pvalues","g","edesign", "groups.vector")
        } 
        names(sig.genes) <- c("independ", colnames(dis))
        summary <- as.data.frame(summary[apply(summary,1,not.all.empty),])
## Vars == groups		
      } else if (vars == "groups") {
        if (is.null(groups.vector)) {
          if(is.null(tstep$groups.vector)) {
            stop ("groups.code is missing")
          } else  {
            groups.vector <- tstep$groups.vector
          }		
        }
	groups.vector2 <- c(groups.vector[1], groups.vector)
        group <- unique(groups.vector2)
        summary <- matrix(" ", ncol = length(group), nrow = nrow(sig.profiles))
        colnames (summary) <- group
        sig.genes <- as.list(group)
        for (i in 1:length(group)) {
          group.sig <- sig.pvalues[,3:(length(groups.vector2)+2)]
          group.sig <- group.sig[,groups.vector2 == group[i]]
          ps <- sig.pvalues[which(apply(group.sig, 1, any)),]
          sigs <- sig.profiles[which(apply(group.sig, 1, any)),]
          coeffs <- coefficients[which(apply(group.sig, 1, any)),]
          ts <- t.score[which(apply(group.sig, 1, any)),]
          if (nrow(sigs) > 0) 
            names.sigs <- rownames(sigs) 
          else names.sigs <- NULL
          summary[,i] <- c(names.sigs, rep(" ", nrow(sig.profiles)-nrow(sigs)))
          sig.genes [[i]] <- list(sigs, coeffs, ts,ps,nrow(ps), edesign, groups.vector)
          names(sig.genes[[i]]) <- c("sig.profiles", "coefficients","t.score", 
                "sig.pvalues", "g", "edesign", "groups.vector")
        }
        names(sig.genes) <- unique(groups.vector2) 
        if (nrow(summary) > 1)  summary <- as.data.frame(summary[apply(summary,1,not.all.empty),])
      } else stop("invalid vars value, must be one of: all, each, groups")
##  Trat. repl.spots
      if(trat.repl.spots == "average") {
        if(vars != "all") {
          for (i in 1:length(sig.genes)) {
            sig.genes[[i]][[1]] <- average.rows(sig.genes[[i]][[1]],index=index, match=match,r=r)
            sig.genes[[i]][[2]] <- average.rows(sig.genes[[i]][[2]],index=index, match=match,r=0)
            sig.genes[[i]][[2]] <- sig.genes[[i]][[2]][is.element (sig.genes[[i]][[2]],sig.genes[[i]][[1]]),]
            sig.genes[[i]][[3]] <- average.rows(sig.genes[[i]][[3]],index=index, match=match,r=r)
            sig.genes[[i]][[3]] <- sig.genes[[i]][[3]][is.element (sig.genes[[i]][[3]],sig.genes[[i]][[1]]),]
          }
        } else {
          sig.genes[[1]] <- average.rows(sig.genes[[1]],index=index, match=match,r=r)
          sig.genes[[2]] <- average.rows(sig.genes[[2]],index=index, match=match,r=r)
          sig.genes[[2]] <- sig.genes[[2]][is.element (sig.genes[[2]],sig.genes[[1]]),]
          sig.genes[[3]] <- average.rows(sig.genes[[3]],index=index, match=match,r=r)
          sig.genes[[2]] <- sig.genes[[3]][is.element (sig.genes[[3]],sig.genes[[1]]),]
        }
      }	 
## Only names
    sig.genes2 <- sig.genes
     if(only.names && vars != "all") {
       for( i in 1 : length(sig.genes)) {   
         if(!is.null(dim(sig.genes[[i]][[1]]))) {  
           sig.genes[[i]][[1]] <- rownames(sig.genes[[i]][[1]])
         }
       }
     }
## Add gene IDs
    if(add.IDs && vars != "all") { 
      for( i in 1 : length(sig.genes)) {     
        if(nrow(sig.genes2[[i]][[1]]) > 1 ) { 
          row.names <- rownames(sig.genes2[[i]][[1]])
          if(trat.repl.spots == "none") {
            ids <- IDs[is.element(IDs[,matchID.col], rownames(sig.genes2[[i]][[1]])),]
          } else { stop("function parameters no compatible (add.IDs, trat.repl.spots)") }
          if(!only.names) {
            sig.genes[[i]][[1]] <- cbind(ids, sig.genes[[i]][[1]])
          } else sig.genes[[i]][[1]] <- ids 
          rownames(sig.genes[[i]][[1]]) <- row.names
        } 
      }
    }
  } else {
     sig.genes <- NULL
     summary <- c("no significant genes")
     print ("no significant genes")
  }
## Output
  output  <- list(sig.genes, summary)
  names(output) <- c("sig.genes", "summary")
  output
}

