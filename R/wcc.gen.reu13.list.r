### For wallace format.

gen.reu13.list <- function(seq.string, aa.list = .CF.GV$amino.acid,
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
  ret <- build.reu13.list(seq.string, aa.list, split.S = split.S)

  ret
} # End of gen.reu13.list().


build.reu13.list <- function(seq.string, aa.list, split.S = TRUE){
  names.seq <- names(seq.string)

  if(split.S){
    synonymous.codon <- .CF.GV$synonymous.codon.split
  } else{
    synonymous.codon <- .CF.GV$synonymous.codon
  }

  ret <- list()
  for(i.gene in 1:length(seq.string)){
    ret[[i.gene]] <- list()
    for(i.aa in 1:length(aa.list)){
      scodon <- synonymous.codon[[aa.list[i.aa]]]
      ret[[i.gene]][[i.aa]] <- list()
      for(i.codon in 1:length(scodon)){
        ret[[i.gene]][[i.aa]][[i.codon]] <-
          which(seq.string[[i.gene]] %in% scodon[i.codon])
      }
      names(ret[[i.gene]][[i.aa]]) <- scodon
    }
    names(ret[[i.gene]]) <- aa.list
  }
  names(ret) <- names.seq

  ret
} # End of build.reu13.list().

