### For wallace format.

gen.y <- function(seq.string, aa.list = .CF.GV$amino.acid,
    split.S = TRUE, drop.X = TRUE, drop.MW = TRUE){
  if(split.S){
    if("S" %in% aa.list){
      if(! "Z" %in% aa.list){
        aa.list <- c(aa.list, "Z")
      }
    } else{
      split.S <- FALSE
    }
  } else{
    if(all(c("S", "Z") %in% aa.list)){
        split.S <- TRUE
    }
  }

  if(drop.X){
    aa.list <- aa.list[aa.list != "X"]
  }

  if(drop.MW){
    aa.list <- aa.list[!(aa.list %in% c("M", "W"))]
  }

  aa.list <- sort(aa.list)
  ret <- build.y(seq.string, aa.list, split.S = split.S)

  ret <- rearrange.y(ret)
  ret
} # End of gen.y().


build.y <- function(seq.string, aa.list, split.S = TRUE){
  names.seq <- names(seq.string)

  if(split.S){
    synonymous.codon <- .CF.GV$synonymous.codon.split
  } else{
    synonymous.codon <- .CF.GV$synonymous.codon
  }

  ret <- list()
  for(i.aa in 1:length(aa.list)){
    scodon <- synonymous.codon[[aa.list[i.aa]]]
    tmp <- lapply(1:length(seq.string),
             function(i.gene){
               tmp.ret <- rep(0L, length(scodon))
               for(i.codon in 1:length(scodon)){
                 tmp.ret[i.codon] <- sum(seq.string[[i.gene]] == scodon[i.codon])
               }
               as.integer(tmp.ret)
             })
    ret[[i.aa]] <- do.call("rbind", tmp)
    colnames(ret[[i.aa]]) <- scodon
    rownames(ret[[i.aa]]) <- names.seq
  }
  names(ret) <- aa.list

  ret
} # End of build.y().

