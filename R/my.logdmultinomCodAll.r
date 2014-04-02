# This compute log posterior on all amino acids.

my.logdmultinomCodAllR <- function(b, phi, y, n, reu13.df = NULL){
  # Returns log posterior of codon draws for all amino acids.
  # For each element, it is a vector of length "# of genes".
  lpclist <- lapply(1:length(y),
               function(i.aa){ # i'th amino acid.
                 .cubfitsEnv$my.logdmultinomCodOne(
                   b[[i.aa]], phi, y[[i.aa]], n[[i.aa]], vec = TRUE,
                   reu13.df.aa = reu13.df[[i.aa]])
               })

  # Return posterior which is the sum of all amino acids.
  rowSums(do.call("cbind", lpclist))
} # End of my.logdmultinomCodAllR().
