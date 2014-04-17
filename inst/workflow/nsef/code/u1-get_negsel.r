### Convert delta t to negative selection.

get.negsel <- function(b.PM, id.slop, aa.names, b.label, b.ci.PM = NULL){
  b.negsel.PM <- b.PM[id.slop]
  if(!is.null(b.ci.PM)){
    b.negsel.ci.PM <- matrix(b.ci.PM[id.slop,], ncol = 2)
  }
  b.negsel.label <- b.label
  for(i.aa in aa.names){
    id.aa <- grep(paste("^", i.aa, "\\.", sep = ""), b.label)

    if(any(b.negsel.PM[id.aa] > 0)){
      tmp <- max(b.negsel.PM[id.aa])
      id.max <- which.max(b.negsel.PM[id.aa])
      b.negsel.PM[id.aa] <- b.negsel.PM[id.aa] - tmp
      b.negsel.PM[id.aa][id.max] <- -tmp

      if(!is.null(b.ci.PM)){
        if(length(id.aa) == 1){
          b.negsel.ci.PM[id.aa,] <- -b.negsel.ci.PM[id.aa,]
        } else{
          tmp.ci <- matrix(b.negsel.ci.PM[id.aa,], ncol = 2)[id.max,]
          b.negsel.ci.PM[id.aa,] <- b.negsel.ci.PM[id.aa,] - tmp
          b.negsel.ci.PM[id.aa,][id.max,] <- -tmp.ci
        }
      }

      # Replace codon name.
      tmp <- .CF.GV$synonymous.codon[[i.aa]]
      if(sum(id.aa) != length(tmp)){
        tmp <- .CF.GV$synonymous.codon.split[[i.aa]]
      }
      b.negsel.label[id.aa][id.max] <-
        paste(i.aa, ".", tmp[length(tmp)], sep = "")
    }
  }

  if(!is.null(b.ci.PM)){
    ret <- list(b.negsel.PM = b.negsel.PM, b.negsel.ci.PM = b.negsel.ci.PM,
                b.negsel.label = b.negsel.label)
  } else{
    ret <- list(b.negsel.PM = b.negsel.PM,
                b.negsel.label = b.negsel.label)
  }

  ret
} # End of get.negsel().
