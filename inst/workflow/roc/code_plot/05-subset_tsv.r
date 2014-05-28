### This script collect the poster means of MCMC runs.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))
source("00-set_env.r")

### Get all cases.
for(i.case in case.names){
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    next
  }
  load(fn.in)
  fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  if(!file.exists(fn.in)){
    next
  }
  load(fn.in)

  ### For log.mu.
  all.names <- names(b.logmu.PM)
  id.var <- grep("log.mu", all.names)

  AA <- gsub("(.)\\.(.*)", "\\1", b.logmu.label)
  CODON <- gsub("(.)\\.(.*)", "\\2", b.logmu.label)

  tmp.PM <- b.logmu.PM[id.var]
  tmp.CI <- b.logmu.ci.PM[id.var,]
  ret <- data.frame(AA = AA, CODON = CODON, Mean = tmp.PM,
                    CI.025 = tmp.CI[, 1], CI.975 = tmp.CI[, 2])
  for(i.aa in unique(AA)){
    tmp <- .CF.GV$synonymous.codon.split[[i.aa]]
    tmp <- data.frame(AA = i.aa, CODON = tmp[length(tmp)],
                      Mean = 0, CI.025 = 0, CI.975 = 0)
    ret <- rbind(ret, tmp)
  }
  order.id <- order(ret$AA, ret$CODON)
  ret <- ret[order.id,] 

  fn.out <- paste(prefix$table, "logmu_", i.case, "_PM.tsv", sep = "")
  write.table(ret, file = fn.out, quote = FALSE, sep = "\t", row.names = FALSE)

  ### For Delta.t.
  all.names <- names(b.negsel.PM)
  id.var <- grep("Delta.t", all.names)

  AA <- gsub("(.)\\.(.*)", "\\1", b.negsel.label)
  CODON <- gsub("(.)\\.(.*)", "\\2", b.negsel.label)

  tmp.PM <- b.negsel.PM[id.var]
  tmp.CI <- b.negsel.ci.PM[id.var,]
  ret <- data.frame(AA = AA, CODON = CODON, Mean = tmp.PM,
                    CI.025 = tmp.CI[, 1], CI.975 = tmp.CI[, 2])
  for(i.aa in unique(AA)){
    tmp <- .CF.GV$synonymous.codon.split[[i.aa]]
    tmp <- data.frame(AA = i.aa, CODON = tmp[length(tmp)],
                      Mean = 0, CI.025 = 0, CI.975 = 0)
    ret <- rbind(ret, tmp)
  }
  order.id <- order(ret$AA, ret$CODON)
  ret <- ret[order.id,] 

  fn.out <- paste(prefix$table, "deltat_", i.case, "_PM.tsv", sep = "")
  write.table(ret, file = fn.out, quote = FALSE, sep = "\t", row.names = FALSE)

  ### For E[Phi].
  ret <- data.frame(ORF = names(phi.PM), Mean = phi.PM, Median = phi.MED,
                    CI.025 = phi.CI[, 1], CI.975 = phi.CI[, 2])

  fn.out <- paste(prefix$table, "phi_", i.case, "_PM.tsv", sep = "")
  write.table(ret, file = fn.out, quote = FALSE, sep = "\t", row.names = FALSE)
}
