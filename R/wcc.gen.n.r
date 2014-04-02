### For wallace format.

gen.n <- function(seq.string, aa.list = .CF.GV$amino.acid,
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
  ret <- build.n(seq.string, aa.list, split.S = split.S)

  ret <- rearrange.n(ret)
  ret
} # End of gen.n().


build.n <- function(seq.string, aa.list, split.S = TRUE){
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
               sum(seq.string[[i.gene]] %in% scodon)
             })
    ret[[i.aa]] <- as.integer(do.call("c", tmp))
    names(ret[[i.aa]]) <- names.seq
  }
  names(ret) <- aa.list

  ret
} # End of build.n().

