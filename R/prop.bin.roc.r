### Wrapper of find.prop.bin.roc().
prop.bin.roc <- function(reu13.df, phi.Obs = NULL, nclass = 20){
  if(is.null(phi.Obs)){
    n.aa <- length(reu13.df)
    tmp.ORF <- lapply(1:n.aa, function(i.aa) unique(reu13.df[[i.aa]]$ORF))
    tmp.ORF <- sort(unique(do.call("c", tmp.ORF)))
    tmp.phi.Obs <- lapply(1:length(tmp.ORF),
                       function(i.orf){
                         for(i.aa in 1:n.aa){
                           id <- which(reu13.df[[i.aa]]$ORF == tmp.ORF[i.orf])
                           if(length(id) > 0){
                             tmp <- reu13.df[[i.aa]]$phi[id[1]]
                             break
                           }
                         }
                         tmp
                       })
    phi.Obs <- do.call("c", tmp.phi.Obs)
    names(phi.Obs) <- tmp.ORF
  }

  ### Rebuild phi.method from phi.Obs or reu13.df
  phi.method <- data.frame(ORF = names(phi.Obs), phi = phi.Obs,
                           stringsAsFactors = FALSE)

  ### Call find.prop.bin.roc().
  phi.bin <- as.vector(quantile(phi.Obs, c(0, seq(0.05, 0.95, length = nclass), 1)))
  ret <- find.prop.bin.roc(reu13.df, phi.method, phi.bin)
  ret
} # End of prop.bin.roc().


### Summarize by amino acid.
### sub.aa.dfs is from REU12.
### phi.method is a data.frame(ORF, phi),
###            where ORF is the key to match codons in sub.aa.dfs.
### phi.bin is a vector used to bin phi.method$phi, usually
###       quantile(phi.method$phi, c(0, seq(0.05, 0.95, length = 19), 1)).

find.prop.bin.roc <- function(sub.aa.dfs, phi.method, phi.bin){
  nclass <- length(phi.bin) - 1
  phi <- phi.method$phi

  ret <- list()
  for(i.aa in 1:length(sub.aa.dfs)){
    tmp.aa <- sub.aa.dfs[[i.aa]]
    u.codon <- as.character(unique(tmp.aa$Codon))

    ### For leftest 2
    i.bin <- 2
    tmp.id <- tmp.aa$ORF %in% phi.method$ORF[phi < phi.bin[i.bin]]
    if(sum(tmp.id) > 0){
      ret.aa <- find.prop(tmp.id, tmp.aa)
      ret.aa$center <- mean(phi.bin[(i.bin - 1):i.bin])
    }

    ### For between 2 to nclass - 1
    for(i.bin in 2:(nclass - 1)){
      tmp.id <- tmp.aa$ORF %in%
                phi.method$ORF[phi > phi.bin[i.bin] & phi < phi.bin[i.bin + 1]]
      if(sum(tmp.id) > 0){
        tmp.prop <- find.prop(tmp.id, tmp.aa)
        tmp.prop$center <- mean(phi.bin[i.bin:(i.bin + 1)])
        ret.aa <- rbind(ret.aa, tmp.prop)
      }
    }

    ### For rightest nclass
    i.bin <- nclass
    tmp.id <- tmp.aa$ORF %in% phi.method$ORF[phi > phi.bin[i.bin]]
    if(sum(tmp.id) > 0){
      tmp.prop <- find.prop(tmp.id, tmp.aa)
      tmp.prop$center <- mean(phi.bin[i.bin:(i.bin + 1)])
      ret.aa <- rbind(ret.aa, tmp.prop)
    }

    ### Bind to all
    ret[[i.aa]] <- ret.aa
  }

  names(ret) <- names(sub.aa.dfs)
  ret
} # End of find.prop.bin.roc().


find.prop <- function(tmp.id, tmp.aa){
  ### tmp.aa is in parent's environment.
  tmp.aa.bin <- tmp.aa[tmp.id,]
  u.orf <- as.character(unique(tmp.aa.bin$ORF))

  ### This is averaged by genes.
  # tmp.table.orf <- NULL
  # for(i.orf in u.orf){
  #   tmp.table <- table(tmp.aa.bin$Codon[tmp.aa.bin$ORF == i.orf])
  #   tmp.table.orf <- rbind(tmp.table.orf, tmp.table / sum(tmp.table))
  # }
  ### This is much fast.
  tmp.table <- table(tmp.aa.bin$ORF, tmp.aa.bin$Codon)
  tmp.table.orf <- tmp.table / rowSums(tmp.table)
  tmp.table.orf <- matrix(tmp.table.orf, nrow = nrow(tmp.table))

  mean.p <- colMeans(tmp.table.orf)
  std.p <- apply(tmp.table.orf, 2, sd)
  stderr.p <- std.p / sum(tmp.id)

  name.codon <- colnames(tmp.table)
  ret <- data.frame(codon = as.character(name.codon),
                    freq.mean = as.double(mean.p),
                    freq.std = as.double(std.p),
                    freq.stderr = as.double(stderr.p),
                    stringsAsFactors = FALSE)
  ret
} # End of find.prop().

